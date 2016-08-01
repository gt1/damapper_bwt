# damapper_bwt
damapper plus BWT based seed searching

Source
------

The damapper_bwt source code is hosted on github:

        https://github.com/gt1/damapper_bwt

It can be obtained using

```
git clone --recursive https://github.com/gt1/damapper_bwt
```

Compilation of damapper_bwt
---------------------------

damapper_bwt needs libmaus2 [https://github.com/gt1/libmaus2].
When libmaus2 is installed in ${LIBMAUSPREFIX} then damapper_bwt can be compiled and
installed in ${HOME}/damapper_bwt using

        - autoreconf -i -f
        - ./configure --with-libmaus2=${LIBMAUSPREFIX} --with-daligner=${DALIGNERLIBPREFIX} \
                --prefix=${HOME}/damapper_bwt
        - make install

