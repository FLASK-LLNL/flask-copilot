FROM python:3.11-bookworm

RUN curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.40.3/install.sh | bash
RUN . $HOME/.nvm/nvm.sh && nvm install 22

WORKDIR /app
COPY requirements.txt /app
COPY flask-app/. /app
COPY charge-backend /app/charge-backend

RUN . $HOME/.nvm/nvm.sh && \
    npm install && \
    npm run build

RUN pip install -r requirements.txt
RUN git clone --recursive https://github.com/FLASK-LLNL/ChARGe.git charge && \
    cd /app/charge && \
    pip install . && \
    charge-install --extras aizynthfinder --extras autogen --extras rdkit --extras chemprice

COPY mock_server.py /app

RUN chmod -R g+rx /app

ENV FLASK_APPDIR=/app/dist

EXPOSE 8001
CMD ["--host", "0.0.0.0", "--port", "8001", "--workers", "8", "mock_server:app"]
ENTRYPOINT ["uvicorn"]

