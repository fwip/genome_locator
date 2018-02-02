FROM python:3.6.4

COPY GRCh38_no_alts.2bit /

COPY requirements.txt /
RUN pip3 install -r requirements.txt
COPY *py /

ENV LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    PIPENV_HIDE_EMOJIS=1

RUN ./build.py GRCh38_no_alts.2bit 10 10

ENV FLASK_APP=run.py

EXPOSE 5000

CMD flask run -h 0.0.0.0 -p 5000
