# Use the official Ubuntu base image
FROM ubuntu:latest

# Set the maintainer label
LABEL maintainer="david.topping@manchester.ac.uk"

# Update the package list and install required packages
RUN apt-get update && apt-get install -y \
    git \
    curl \
    bzip2 \
    && rm -rf /var/lib/apt/lists/*

# Download the Anaconda installer script
RUN curl -O https://repo.anaconda.com/archive/Anaconda3-2023.03-Linux-x86_64.sh

# Install Anaconda
RUN bash Anaconda3-2023.03-Linux-x86_64.sh -b -p /opt/anaconda

# Set the PATH environment variable to include the Anaconda installation
ENV PATH=/opt/anaconda/bin:$PATH

# Clean up
RUN rm Anaconda3-2023.03-Linux-x86_64.sh

# add the repo host to know_hosts, only run once
RUN mkdir -p -m 0700 ~/.ssh && ssh-keyscan github.com >> ~/.ssh/known_hosts

# Clone the GitHub repository
RUN --mount=type=ssh,id=github \ 
git clone git@github.com:loftytopping/AeroSysProp.git /opt/AeroSysProp

# Set the working directory
WORKDIR /opt/AeroSysProp

# Create the conda environment from the environment file
RUN conda env create -f /opt/AeroSysProp/environment.yml

# Activate the conda environment and ensure it's available in the Docker container
RUN echo "source activate $(head -1 /opt/AeroSysProp/environment.yml | cut -d ' ' -f 2)" > ~/.bashrc

# docker build --ssh github=~/.ssh/id_rsa_github -t your-image-name .

# (Optional) Run any setup scripts or commands required by your repository
# RUN ./setup.sh

# Expose any ports your application needs
# EXPOSE 8080

# Define the default command to run your application
CMD ["bash"]


