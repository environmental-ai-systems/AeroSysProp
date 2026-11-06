"""Compare the two guidance levels in a marked campaign.

    python -m evaluation.analyse evaluation/results/2026-10-02-full

Reads marks.csv and writes analysis/ beside it: one CSV per table and
summary.md. Only the automatic marks are used; blind-marked criteria are not.

The comparison is paired. Every model meets every question under both levels,
so the effect for one model on one question is the mean over repeats at `full`
minus the mean at `bare`. Questions, not transcripts, are the unit of
resampling: repeats of one question are not independent, and the question set
is a sample of what users might ask. Intervals are 95% percentile bootstraps
over questions (and, for the pooled effect, over models as well). A pair with a
harness error on either side is dropped, not counted as a fail.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

N_BOOT = 4000
SEED = 20261023

# Criteria that make up the overall automatic score. value_key_all repeats
# value_key at a finer grain, so only one of them counts.
SCORE_CRITERIA = ['tool_called', 'value_key', 'method_ok', 'method_named', 'value_with_units',
                  'version_named', 'caveat_passed', 'no_false_warning', 'failure_reported',
                  'batch_tool_used', 'run_folder_cited', 'coverage_stated', 'basis_named',
                  'status_reported']

# Outcome name -> (description, which rows it applies to).
OUTCOMES = {
    'score': 'Mean of the applicable automatic criteria',
    'tool_called': 'At least one tool call the server accepted',
    'value_key': 'Fraction of key values found in the answer',
    'caveat_passed': 'Diagnostic caveat passed on to the user',
    'no_false_warning': 'No warning given where none was due',
    'batch_steps': 'Batch tool used, run folder cited, coverage stated (mean)',
    'invention': 'Answer has numbers no tool returned (lower is better)',
}

# Parameters in billions, and active parameters for mixture-of-experts models.
ACTIVE = {'qwen3:30b-a3b': 3.3, 'gpt-oss:20b': 3.6}


def size_b(model: str) -> float:
    m = re.search(r':(\d+(?:\.\d+)?)b', model)
    return float(m.group(1)) if m else float('nan')


JUDGE_COUNTS = ['j_numbers_checked', 'j_invented', 'j_over_rounded', 'j_unclear',
                'j_has_numbers', 'j_judged']


def add_judge(m: pd.DataFrame, judge_marks: pd.DataFrame) -> tuple:
    """Join judge marks (evaluation.judge merge) onto the automatic ones.

    judge_score is the share of judged criteria met. judge_invented is 1 where the judge
    called any number invented and 0 where it found none; a completed transcript with no
    number lacking a successful tool result behind it counts as 0. judge_over_rounded is the
    same for values rounded more coarsely than the marking tolerance, and judge_unsupported
    for either, the headline: no question asks for rounded values. Transcripts whose numbers
    are still waiting for the judge are left out."""
    keys = ['model', 'guidance', 'question_id', 'repeat']
    jm = judge_marks[judge_marks['j_judged'] == 1]
    crit = [c for c in jm.columns if c.startswith('j_') and c not in JUDGE_COUNTS]
    m = m.merge(jm[keys + crit + [c for c in JUDGE_COUNTS if c in jm]], on=keys, how='left')
    m['judge_score'] = m[crit].mean(axis=1, skipna=True) if crit else np.nan
    pending = judge_marks[(judge_marks['j_judged'] == 0) & (judge_marks['j_has_numbers'] == 1)]
    waiting = m.set_index(keys).index.isin(pending.set_index(keys).index)
    with_numbers = m['j_has_numbers'].fillna(0) == 1
    if 'j_over_rounded' not in m:  # judged before over_rounded existed
        m['j_over_rounded'] = np.nan

    def flag(counts):
        out = pd.Series(0.0, index=m.index).where(m['status'] == 'complete')
        out[with_numbers] = (counts[with_numbers].fillna(0) > 0).astype(float)
        out[waiting] = np.nan
        return out

    m['judge_invented'] = flag(m['j_invented'])
    m['judge_over_rounded'] = flag(m['j_over_rounded'])
    m['judge_unsupported'] = flag(m['j_invented'].fillna(0) + m['j_over_rounded'].fillna(0))
    extra = {'judge_score': 'Judge: share of judged criteria met',
             'judge_unsupported': 'Judge: answer contains an invented or over-rounded value '
                                  '(lower is better)',
             'judge_invented': 'Judge: answer contains an invented value (lower is better)',
             'judge_over_rounded': 'Judge: answer rounds a value beyond the tolerance '
                                   '(lower is better)'}
    extra.update({c: f'Judge: {c[2:]}' for c in crit})
    return m, extra


def outcomes(marks: pd.DataFrame) -> pd.DataFrame:
    m = marks.copy()
    m['score'] = m[SCORE_CRITERIA].mean(axis=1, skipna=True)
    batch = ['batch_tool_used', 'run_folder_cited', 'coverage_stated']
    m['batch_steps'] = m[batch].mean(axis=1, skipna=True).where(m['kind'] == 'batch')
    m['invention'] = (m['unsupported_numbers'] > 0).astype(float).where(
        m['unsupported_numbers'].notna())
    return m


def paired(m: pd.DataFrame, outcome: str) -> pd.DataFrame:
    """One row per model and question: mean over repeats at each level, and the difference."""
    d = m[m['status'] == 'complete'].dropna(subset=[outcome])
    cell = d.groupby(['model', 'question_id', 'kind', 'guidance'])[outcome].mean().unstack()
    if not {'full', 'bare'} <= set(cell.columns):  # nothing to pair yet, e.g. a partial judge run
        return pd.DataFrame(columns=['model', 'question_id', 'kind', 'bare', 'full', 'diff'])
    cell = cell.dropna(subset=['full', 'bare']).reset_index()
    cell['diff'] = cell['full'] - cell['bare']
    return cell


def _boot_questions(cell: pd.DataFrame, rng, col='diff') -> np.ndarray:
    """Bootstrap over questions; each resample averages per question, then over questions."""
    per_q = cell.groupby('question_id')[col].mean().to_numpy()
    idx = rng.integers(0, len(per_q), (N_BOOT, len(per_q)))
    return per_q[idx].mean(axis=1)


def _boot_two_way(cell: pd.DataFrame, rng) -> np.ndarray:
    """Bootstrap over models and questions together, for an effect pooled across models."""
    table = cell.pivot_table(index='model', columns='question_id', values='diff')
    a = table.to_numpy()
    nm, nq = a.shape
    out = np.empty(N_BOOT)
    for b in range(N_BOOT):
        s = a[np.ix_(rng.integers(0, nm, nm), rng.integers(0, nq, nq))]
        out[b] = np.nanmean(np.nanmean(s, axis=1))
    return out


def _ci(samples) -> tuple:
    if len(samples) == 0:  # e.g. a correlation over two models, or sizes unknown
        return float('nan'), float('nan')
    lo, hi = np.percentile(samples, [2.5, 97.5])
    return round(float(lo), 3), round(float(hi), 3)


def _p_sign(samples) -> float:
    """Two-sided bootstrap p for a mean of zero: twice the smaller tail."""
    s = np.asarray(samples)
    return round(float(min(1.0, 2 * min((s <= 0).mean(), (s >= 0).mean()))), 4)


def pooled_effects(m: pd.DataFrame, rng, outcomes_=None) -> pd.DataFrame:
    rows = []
    for name, desc in (outcomes_ or OUTCOMES).items():
        cell = paired(m, name)
        if cell.empty:
            continue
        per_model = cell.groupby('model')[['bare', 'full', 'diff']].mean()
        boot = _boot_two_way(cell, rng)
        rows.append({'outcome': name, 'description': desc,
                     'models': cell['model'].nunique(), 'questions': cell['question_id'].nunique(),
                     'bare': round(per_model['bare'].mean(), 3),
                     'full': round(per_model['full'].mean(), 3),
                     'diff': round(per_model['diff'].mean(), 3),
                     'ci_low': _ci(boot)[0], 'ci_high': _ci(boot)[1], 'p': _p_sign(boot),
                     'models_better': int((per_model['diff'] > 0).sum()),
                     'models_worse': int((per_model['diff'] < 0).sum())})
    return pd.DataFrame(rows)


def per_model(m: pd.DataFrame, outcome: str, rng) -> pd.DataFrame:
    cell = paired(m, outcome)
    rows = []
    for model, c in cell.groupby('model'):
        boot = _boot_questions(c, rng)
        q = c.groupby('question_id')['diff'].mean()
        rows.append({'model': model, 'params_b': size_b(model),
                     'active_b': ACTIVE.get(model, size_b(model)),
                     'bare': round(c['bare'].mean(), 3), 'full': round(c['full'].mean(), 3),
                     'diff': round(c['diff'].mean(), 3), 'ci_low': _ci(boot)[0],
                     'ci_high': _ci(boot)[1], 'p': _p_sign(boot),
                     'questions_better': int((q > 0).sum()), 'questions_worse': int((q < 0).sum())})
    out = pd.DataFrame(rows)
    return out.sort_values('diff', ascending=False) if rows else out


def per_kind(m: pd.DataFrame, outcome: str, rng) -> pd.DataFrame:
    cell = paired(m, outcome)
    rows = []
    for kind, c in cell.groupby('kind'):
        boot = _boot_two_way(c, rng) if c['question_id'].nunique() > 1 else [np.nan]
        pm = c.groupby('model')[['bare', 'full', 'diff']].mean()
        rows.append({'kind': kind, 'questions': c['question_id'].nunique(),
                     'bare': round(pm['bare'].mean(), 3), 'full': round(pm['full'].mean(), 3),
                     'diff': round(pm['diff'].mean(), 3),
                     'ci_low': _ci(boot)[0] if len(boot) > 1 else np.nan,
                     'ci_high': _ci(boot)[1] if len(boot) > 1 else np.nan})
    return pd.DataFrame(rows).sort_values('diff', ascending=False)


def per_question(m: pd.DataFrame, outcome: str) -> pd.DataFrame:
    cell = paired(m, outcome)
    q = cell.groupby(['question_id', 'kind']).agg(
        bare=('bare', 'mean'), full=('full', 'mean'), diff=('diff', 'mean'),
        models_better=('diff', lambda s: int((s > 0).sum())),
        models_worse=('diff', lambda s: int((s < 0).sum()))).round(3).reset_index()
    return q.sort_values('diff', ascending=False)


def per_kind_model(m: pd.DataFrame, outcome: str) -> pd.DataFrame:
    cell = paired(m, outcome)
    return cell.pivot_table(index='model', columns='kind', values='diff', aggfunc='mean').round(3)


def repeat_agreement(m: pd.DataFrame, outcome: str) -> pd.DataFrame:
    """How often the three repeats agree, per model and level: a guide to how noisy a cell is."""
    d = m[m['status'] == 'complete'].dropna(subset=[outcome])
    g = d.groupby(['model', 'guidance', 'question_id'])[outcome]
    spread = (g.max() - g.min()).reset_index(name='spread')
    return (spread.groupby(['model', 'guidance'])['spread']
            .agg(lambda s: round(float((s == 0).mean()), 3)).unstack())


def _rank_corr(x, y) -> float:
    rx, ry = pd.Series(x).rank(), pd.Series(y).rank()
    return round(float(np.corrcoef(rx, ry)[0, 1]), 3)


def who_benefits(model_table: pd.DataFrame, rng) -> dict:
    """Spearman correlations of the effect with size and with the bare score, with bootstrap
    intervals over models. A negative correlation with the bare score is partly regression to
    the mean and the ceiling, so it is reported beside size rather than instead of it."""
    t = model_table.reset_index(drop=True)
    out = {}
    for col in ('params_b', 'active_b', 'bare'):
        boots = []
        for _ in range(N_BOOT):
            s = t.iloc[rng.integers(0, len(t), len(t))]
            if s[col].nunique() > 1 and s['diff'].nunique() > 1:
                boots.append(_rank_corr(s[col], s['diff']))
        out[col] = (_rank_corr(t[col], t['diff']), *_ci(boots))
    return out


def data_quality(answers: pd.DataFrame) -> pd.DataFrame:
    """Per model and level: empty final answers, thinking left in the answer (which the marker
    reads as if it were the answer), and question cells whose repeats are byte-identical
    (repeats share a seed, so they are not independent draws)."""
    a = answers[answers['status'] == 'complete'].copy()
    a['fa'] = a['final_answer'].fillna('')
    a['empty'] = a['fa'].str.strip() == ''
    a['think_in_answer'] = a['fa'].str.contains('<think>', regex=False)
    flags = a.groupby(['model', 'guidance'])[['empty', 'think_in_answer']].mean()
    same = (a.groupby(['model', 'guidance', 'question_id'])['fa'].nunique() == 1)
    flags['repeats_identical'] = same.groupby(['model', 'guidance']).mean()
    return flags.unstack().round(3)


def _md(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    lines = ['| ' + ' | '.join(map(str, cols)) + ' |', '|' + '---|' * len(cols)]
    for _, r in df.iterrows():
        lines.append('| ' + ' | '.join('' if pd.isna(v) else str(v) for v in r) + ' |')
    return '\n'.join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('campaign', type=Path, help='folder holding marks.csv')
    args = ap.parse_args(argv)
    marks_path = args.campaign / 'marks.csv'
    if not marks_path.exists():
        sys.exit(f'{marks_path} not found: run evaluation.mark first')
    out = args.campaign / 'analysis'
    out.mkdir(exist_ok=True)
    rng = np.random.default_rng(SEED)

    marks = pd.read_csv(marks_path)
    m = outcomes(marks)
    errors = marks[marks['status'] != 'complete']

    judge_path = args.campaign / 'judge' / 'judge_marks.csv'
    judged = None
    if judge_path.exists():
        m, judge_outcomes = add_judge(m, pd.read_csv(judge_path))
        judged = pooled_effects(m, rng, judge_outcomes)
        judged.to_csv(out / 'judge_effects.csv', index=False)
        per_model(m, 'judge_score', rng).to_csv(out / 'per_model_judge_score.csv', index=False)
        per_question(m, 'judge_score').to_csv(out / 'per_question_judge_score.csv', index=False)

    pooled = pooled_effects(m, rng)
    models_score = per_model(m, 'score', rng)
    models_value = per_model(m, 'value_key', rng)
    kinds = per_kind(m, 'score', rng)
    questions = per_question(m, 'score')
    kind_model = per_kind_model(m, 'score')
    agree = repeat_agreement(m, 'score')
    benefit = who_benefits(models_score, rng)
    stops = (m.groupby(['guidance', 'stop_reason']).size().unstack(fill_value=0))
    calls = m[m['status'] == 'complete'].groupby('guidance')[['n_calls', 'tool_errors']].mean().round(2)
    answers_path = args.campaign / 'answers.csv'
    quality = data_quality(pd.read_csv(answers_path)) if answers_path.exists() else None
    if quality is not None:
        quality.columns = [f'{c}_{g}' for c, g in quality.columns]
        quality.to_csv(out / 'data_quality.csv')

    pooled.to_csv(out / 'pooled_effects.csv', index=False)
    models_score.to_csv(out / 'per_model_score.csv', index=False)
    models_value.to_csv(out / 'per_model_value_key.csv', index=False)
    kinds.to_csv(out / 'per_kind_score.csv', index=False)
    questions.to_csv(out / 'per_question_score.csv', index=False)
    kind_model.to_csv(out / 'per_kind_model_score_diff.csv')
    agree.to_csv(out / 'repeat_agreement_score.csv')
    errors[['model', 'guidance', 'question_id', 'repeat']].to_csv(out / 'harness_errors.csv',
                                                                   index=False)

    benefit_rows = pd.DataFrame([{'against': k, 'spearman': v[0], 'ci_low': v[1], 'ci_high': v[2]}
                                 for k, v in benefit.items()])
    benefit_rows.to_csv(out / 'who_benefits.csv', index=False)

    summary = [
        f'# Guidance effect: {args.campaign.name}', '',
        f'{len(marks)} transcripts, {marks["model"].nunique()} models, '
        f'{marks["question_id"].nunique()} questions. Automatic marks only; '
        f'{len(errors)} harness errors dropped with their pairs (harness_errors.csv).',
        'Effects are full minus bare. Intervals are 95% bootstraps over models and questions.',
        '', '## Pooled effects', '', _md(pooled), '',
        '## Per model (overall automatic score)', '', _md(models_score), '',
        '## Per model (key values found)', '', _md(models_value), '',
        '## Per kind of question (overall score)', '', _md(kinds), '',
        '## Who benefits (Spearman with the score effect, over models)', '', _md(benefit_rows), '',
        '## Questions most and least helped (overall score)', '',
        _md(pd.concat([questions.head(8), questions.tail(8)])), '',
        '## Stop reasons', '', _md(stops.reset_index()), '',
        '## Tool calls per conversation', '', _md(calls.reset_index()), '',
        '## Repeats that agree exactly (overall score)', '', _md(agree.reset_index()), '',
    ]
    if judged is not None:
        summary += ['## Judge marks (evaluation.judge)', '',
                    f'{int(m["judge_score"].notna().sum())} transcripts with judged criteria.', '',
                    _md(judged), '']
    if quality is not None:
        summary += ['## Data quality (share of transcripts; repeats_identical is per question)',
                    '', _md(quality.reset_index()), '']
    (out / 'summary.md').write_text('\n'.join(summary))
    print(f'Wrote {out}/summary.md and {len(list(out.glob("*.csv")))} tables')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
