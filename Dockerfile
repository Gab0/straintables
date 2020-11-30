
# This image sets a environment able to run Python3 with numpy, scipy, pandas,
# pygobject and sklearn. 
#
# It serves as the foundtation on CI tests for straintables.
#
# Image available at dockerhub: gab0/python3-gtk3-scipy:v1



FROM python:3.8.6-alpine3.12

ENV LANG C.UTF-8
RUN echo "http://dl-cdn.alpinelinux.org/alpine/latest-stable/main" > /etc/apk/repositories
RUN echo "http://dl-cdn.alpinelinux.org/alpine/latest-stable/community" >> /etc/apk/repositories
RUN echo "http://dl-8.alpinelinux.org/alpine/edge/community" >> /etc/apk/repositories
RUN apk add --no-cache --allow-untrusted --repository http://dl-3.alpinelinux.org/alpine/edge/testing hdf5 hdf5-dev

# SYSTEM REQUIREMENTS;
RUN apk add --no-cache \
  xauth \
  xvfb \
  xvfb-run \
  gcc \
  gfortran \
  build-base \
  wget \
  freetype-dev \
  libpng-dev \
  openblas-dev \
  glib \
  jpeg-dev \
  zlib-dev \
  jpeg-dev \
  git




# PYTHON REQUIREMENTS;
RUN pip install numpy scipy cython sklearn pandas

# DOWNLOAD AND SETUP CLUSTAL OMEGA;
RUN wget http://www.clustal.org/omega/clustalo-1.2.4-Ubuntu-x86_64 -O clustalo
RUN chmod +x clustalo
RUN mv clustalo /bin/clustalo



COPY . /straintables

# INSTALL STRAINTABLES MODULE & BINARIES;
RUN pip install /straintables
RUN pip install --no-binary :all: ./straintables
RUN stgenomepline --help

ENTRYPOINT ["/bin/sh"]

