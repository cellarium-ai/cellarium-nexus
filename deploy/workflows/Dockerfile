FROM python:3.10-slim-bullseye

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /app

# system deps
RUN apt-get update && apt-get install -y gcc python3-dev libpq-dev \
 && apt-get clean && rm -rf /var/lib/apt/lists/*

# copy code and install it
COPY . /app

RUN pip install --no-cache-dir --upgrade pip setuptools wheel \
 && pip install --no-cache-dir -r requirements/pipelines.txt \
 && pip install --no-cache-dir -e /app

ENV PYTHONPATH=/app:/app/cellarium:/app/cellarium/nexus
