FROM python:3.6.4

RUN useradd -m -U flask
USER flask:flask
WORKDIR /home/flask/

COPY GRCh38_no_alts.2bit .
COPY requirements.txt .
RUN pip3 install --user -r requirements.txt
COPY build_genome_index.py search.py lookup_hash.py ./

# RUN ./build_genome_index.py GRCh38_no_alts.2bit 10 10
COPY GRCh38_no_alts.2bit.M10.Q10.index.hdf5 .

COPY run.py .
ENV FLASK_APP=run.py

EXPOSE 5000

CMD ./run.py
