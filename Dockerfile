# Use an official Ubuntu base image
FROM ubuntu:20.04

# Set non-interactive frontend for apt to avoid prompts during install
ENV DEBIAN_FRONTEND=noninteractive

# Install g++, make, and other essential tools
RUN apt-get update && \
    apt-get install -y \
        build-essential \
        libeigen3-dev \
        g++ \
        cmake \
        git \
        wget \
        vim \
        && apt-get clean && rm -rf /var/lib/apt/lists/*

# Create a working directory
WORKDIR /app/genevol

# Copy your C++ source files into the container
COPY . /app/genevol

# Install dependencies
RUN /app/genevol/script/install_sources.sh --workdir /app

# Build program
RUN /app/genevol/script/build.sh --workdir /app

# Default command to run when starting the container
# (Modify as needed, e.g., "./myapp")
CMD ["/bin/bash"]
