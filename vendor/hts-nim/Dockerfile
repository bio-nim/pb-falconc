FROM alpine:3.11.5
#FROM alpine:20190925

ENV CFLAGS="-fPIC -O3"

RUN apk add wget git xz bzip2-static musl m4 autoconf tar xz-dev bzip2-dev build-base libpthread-stubs libzip-dev gfortran \
	    openssl-libs-static openblas-static pcre-dev curl llvm-dev curl-static bash

RUN mkdir -p /usr/local/include && \
    git clone --depth 1 https://github.com/ebiggers/libdeflate.git && \
    cd libdeflate && make -j4 CFLAGS="-fPIC -O3" install && \
    cd .. && rm -rf libdeflate && \
    git clone https://github.com/cloudflare/zlib cloudflare-zlib && \
    cd cloudflare-zlib && ./configure && make install && \
    cd .. && \
    rm -rf cloudflare-zlib

RUN cd / && \
    git clone -b v1.2.6 git://github.com/nim-lang/nim nim && \
    cd nim &&  \
    sh ./build_all.sh && \
    rm -rf csources && \
    echo 'PATH=/nim/bin:$PATH' >> ~/.bashrc && \
    echo 'PATH=/nim/bin:$PATH' >> ~/.bash_profile && \
    echo 'PATH=/nim/bin:$PATH' >> /etc/environment 

RUN apk add cmake openssl-dev && \
	wget https://libzip.org/download/libzip-1.6.1.tar.gz && \
	tar xzvf libzip-1.6.1.tar.gz && \
	cd libzip-1.6.1 && \
	mkdir build && cd build && \
	cmake -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=/usr/local/ ../ && \
	make -j4 CFLAGS="-fPIC -O3" install && \
	cd ../../ && rm -rf libzip-1.6.1*


ENV PATH=:/root/.nimble/bin:/nim/bin/:$PATH	

RUN \
    git clone https://github.com/samtools/htslib && \
    cd htslib && git checkout 1.11 && autoheader && autoconf && \
    ./configure --disable-s3 --disable-libcurl --with-libdeflate && \
    make -j4 CFLAGS="-fPIC -O3" install && \
    cd ../ && \
    git clone https://github.com/samtools/bcftools && \
    cd bcftools && git checkout 1.10.2 && autoheader && autoconf && \
    ./configure --disable-s3 --disable-libcurl --with-libdeflate && \
    make -j4 CFLAGS="-fPIC -O3" install && \
    cd ../ && rm -rf htslib bcftools


RUN sh -c 'curl --proto "=https" --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -q -y'
#	    && apk add clang-libs

ENV HTSLIB=system
ENV PATH=$PATH:~/.cargo/bin/

#COPY docker/d4.patch /tmp/

#&& git apply < /tmp/d4.patch \
RUN ~/.cargo/bin/rustup target add x86_64-unknown-linux-musl \
	&& git clone https://github.com/38/d4-format \
	&& cd d4-format \
        && ln -s /usr/bin/gcc /usr/bin/musl-gcc \
	&& ~/.cargo/bin/cargo build --all --target x86_64-unknown-linux-musl --release

RUN install -m 644 d4-format/target/x86_64-unknown-linux-musl/release/libd4binding.a /usr/lib && \
	install -m 644 d4-format/d4binding/include/d4.h /usr/include


ADD . /src/
RUN cat /src/docker/docker.nim.cfg >> /nim/config/nim.cfg && \
    source ~/.bashrc && cd /src/ && nimble install -y && \
    nimble install -y c2nim docopt && \
    nimble install -y websocket@#head && \
    nim c -o:/usr/local/bin/nsb /src/docker/nsb.nim && \
    rm -rf /src/

