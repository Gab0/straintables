
# This image sets a environment able to run Python3 with numpy, scipy, pandas,
# pygobject and sklearn. 
#
# It serves as the foundtation on CI tests for straintables.
#
# Image available at dockerhub: gab0/python3-gtk3-scipy:v1



FROM python:3.7.3-alpine

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
  cairo \
  cairo-dev \
  jpeg-dev \
  zlib-dev \
  gobject-introspection-dev \
  py3-cairo \
  py-cairo-dev \
  jpeg-dev \
  git


RUN apk add --no-cache \
  gtk+3.0

# PYTHON REQUIREMENTS;
RUN pip install numpy scipy cython sklearn pandas pygobject


ENTRYPOINT ["/bin/bash"]
#COPY . /app

#RUN pip install /app
#RUN pip install --no-binary :all: ./app
#RUN lmpline --help

