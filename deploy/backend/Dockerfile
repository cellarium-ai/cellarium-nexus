FROM python:3.12-slim-bookworm

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    DJANGO_SETTINGS_MODULE=cellarium.nexus.backend.application.settings

WORKDIR /app

# system deps
RUN apt-get update && apt-get install -y gcc python3-dev libpq-dev git \
 && apt-get clean && rm -rf /var/lib/apt/lists/*

# copy code and install it
COPY . /app

# Make entrypoint script executable in its original path
RUN chmod +x deploy/backend/docker-entrypoint.sh

RUN pip install --no-cache-dir --upgrade pip setuptools wheel \
 && pip install --no-cache-dir -r deploy/requirements/backend.lock.txt \
 && pip install --no-cache-dir . \
 && pip install --no-cache-dir gunicorn

EXPOSE 8080

ENV PYTHONPATH=/app:/app/cellarium:/app/cellarium/nexus

ENTRYPOINT ["deploy/backend/docker-entrypoint.sh"]
