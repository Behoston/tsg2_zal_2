FROM python:3.6

WORKDIR /app
COPY requirements.txt .
RUN pip install pip-tools && pip-sync
COPY . .

CMD ["python", "assembly.py", "--sample"]