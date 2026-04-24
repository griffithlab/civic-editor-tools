FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# Install system packages
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-venv \
    curl \
    wget \
    bash \
    less \
    vim \
    && rm -rf /var/lib/apt/lists/*

# Install NCBI Entrez Direct (efetch, esearch, etc.)
RUN curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh | bash
ENV PATH="/root/edirect:${PATH}"

# Install Python dependencies
COPY requirements.txt /tmp/requirements.txt
RUN pip3 install --break-system-packages -r /tmp/requirements.txt

CMD ["/bin/bash"]
