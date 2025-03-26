FROM genevol_base:latest

# Create a working directory
WORKDIR /app/genevol

# Copy your C++ source files into the container
COPY . /app/genevol

# Build program
RUN /app/genevol/script/build.sh --workdir /app

# Default command to run when starting the container
# (Modify as needed, e.g., "./myapp")
CMD ["/bin/bash"]
