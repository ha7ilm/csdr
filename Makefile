# This software is part of libcsdr, a set of simple DSP routines for 
# Software Defined Radio.
#
# Copyright (c) 2014, Andras Retzler <randras@sdr.hu>
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in the
#      documentation and/or other materials provided with the distribution.
#    * Neither the name of the copyright holder nor the
#      names of its contributors may be used to endorse or promote products
#      derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL ANDRAS RETZLER BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

LIBSOURCES =  fft_fftw.c libcsdr_wrapper.c 
#SOURCES = csdr.c $(LIBSOURCES)
cpufeature = $(if $(findstring $(1),$(shell cat /proc/cpuinfo)),$(2))
PARAMS_SSE = $(call cpufeature,sse,-msse) $(call cpufeature,sse2,-msse2) $(call cpufeature,sse3,-msse3) $(call cpufeature,sse4a,-msse4a) $(call cpufeature,sse4_1,-msse4.1) $(call cpufeature,sse4_2,-msse4.2 -msse4) -mfpmath=sse 
PARAMS_NEON = -mfloat-abi=hard -march=armv7-a -mtune=cortex-a8 -mfpu=neon -mvectorize-with-neon-quad -funsafe-math-optimizations -Wformat=0 -DNEON_OPTS
#tnx Jan Szumiec for the Raspberry Pi support
PARAMS_RASPI = -mfloat-abi=hard -mcpu=arm1176jzf-s -mfpu=vfp -funsafe-math-optimizations -Wformat=0
PARAMS_ARM = $(if $(call cpufeature,BCM2708,dummy-text),$(PARAMS_RASPI),$(PARAMS_NEON))
PARAMS_SIMD = $(if $(call cpufeature,sse,dummy-text),$(PARAMS_SSE),$(PARAMS_ARM))
PARAMS_LOOPVECT = -O3 -ffast-math -fdump-tree-vect-details -dumpbase dumpvect
PARAMS_LIBS = -g -lm -lrt -lfftw3f -DUSE_FFTW -DLIBCSDR_GPL -DUSE_IMA_ADPCM
PARAMS_SO = -fpic  
PARAMS_MISC = -Wno-unused-result
#DEBUG_ON = 0 #debug is always on by now (anyway it could be compiled with `make DEBUG_ON=1`)
#PARAMS_DEBUG = $(if $(DEBUG_ON),-g,)
FFTW_PACKAGE = fftw-3.3.3
PREFIX ?= /usr
PYVERSION ?= 2.7
SOVERSION = 0.15

.PHONY: clean-vect clean-python clean
all: csdr nmux
libcsdr.so: fft_fftw.c fft_rpi.c libcsdr_wrapper.c libcsdr.c libcsdr_gpl.c fastddc.c fastddc.h  fft_fftw.h  fft_rpi.h  ima_adpcm.h  libcsdr_gpl.h  libcsdr.h  predefined.h
	@echo NOTE: you may have to manually edit Makefile to optimize for your CPU \(especially if you compile on ARM, please edit PARAMS_NEON\).
	@echo Auto-detected optimization parameters: $(PARAMS_SIMD)
	@echo
	rm -f dumpvect*.vect
	gcc -std=gnu99 $(PARAMS_LOOPVECT) $(PARAMS_SIMD) $(LIBSOURCES) $(PARAMS_LIBS) $(PARAMS_MISC) -fpic -shared -Wl,-soname,libcsdr.so.$(SOVERSION) -o libcsdr.so.$(SOVERSION)
	@ln -fs libcsdr.so.$(SOVERSION) libcsdr.so
	-./parsevect dumpvect*.vect
csdr: csdr.c libcsdr.so
	gcc -std=gnu99 $(PARAMS_LOOPVECT) $(PARAMS_SIMD) csdr.c $(PARAMS_LIBS) -L. -lcsdr $(PARAMS_MISC) -o csdr
ddcd: ddcd.cpp libcsdr.so ddcd.h
	g++ $(PARAMS_LOOPVECT) $(PARAMS_SIMD) ddcd.cpp $(PARAMS_LIBS) -L. -lcsdr -lpthread $(PARAMS_MISC) -o ddcd
nmux: nmux.cpp libcsdr.so nmux.h tsmpool.cpp tsmpool.h
	g++ $(PARAMS_LOOPVECT) $(PARAMS_SIMD) nmux.cpp tsmpool.cpp $(PARAMS_LIBS) -L. -lcsdr -lpthread $(PARAMS_MISC) -o nmux
arm-cross: clean-vect
	#note: this doesn't work since having added FFTW
	arm-linux-gnueabihf-gcc -std=gnu99 -O3 -fshort-double -ffast-math -dumpbase dumpvect-arm -fdump-tree-vect-details -mfloat-abi=softfp -march=armv7-a -mtune=cortex-a9 -mfpu=neon -mvectorize-with-neon-quad -Wno-unused-result -Wformat=0 $(SOURCES) -lm -o ./csdr
clean-vect:
	rm -f dumpvect*.vect
clean: clean-vect clean-python
	rm -f libcsdr.so.$(SOVERSION) csdr ddcd nmux *.o *.so
install: all 
	install -m 0755 libcsdr.so.$(SOVERSION) $(PREFIX)/lib
	install -m 0755 csdr $(PREFIX)/bin
	install -m 0755 csdr-fm $(PREFIX)/bin
	install -m 0755 nmux $(PREFIX)/bin
	#-install -m 0755 ddcd $(PREFIX)/bin
	@ldconfig || echo please run ldconfig
uninstall:
	rm $(PREFIX)/lib/libcsdr.so.$(SOVERSION) $(PREFIX)/bin/csdr $(PREFIX)/bin/csdr-fm
	ldconfig
disasm:
	objdump -S libcsdr.so.$(SOVERSION) > libcsdr.disasm
emcc-clean:
	-rm sdr.js/sdr.js
	-rm sdr.js/sdrjs-compiled.js
	-rm -rf sdr.js/$(FFTW_PACKAGE)
emcc-get-deps:
	echo "getting and compiling fftw3 with emscripten..."
	cd sdr.js; \
	wget http://fftw.org/$(FFTW_PACKAGE).tar.gz; \
	tar -xvf $(FFTW_PACKAGE).tar.gz; \
	rm $(FFTW_PACKAGE).tar.gz; \
	cd $(FFTW_PACKAGE); \
	emconfigure ./configure --enable-float --disable-fortran --prefix=`pwd`/emscripten-install --libdir=`pwd`/emscripten-lib; \
	emmake make; \
	emmake make install
emcc:
	emcc -O3 -Isdr.js/$(FFTW_PACKAGE)/api -Lsdr.js/$(FFTW_PACKAGE)/emscripten-lib -o sdr.js/sdrjs-compiled.js fft_fftw.c libcsdr_wrapper.c -s TOTAL_MEMORY=67108864 -DLIBCSDR_GPL -DUSE_IMA_ADPCM -DUSE_FFTW -lfftw3f -s EXPORTED_FUNCTIONS="`python sdr.js/exported_functions.py`"
	cat sdr.js/sdrjs-header.js sdr.js/sdrjs-compiled.js sdr.js/sdrjs-footer.js > sdr.js/sdr.js
emcc-beautify:
	bash -c 'type js-beautify >/dev/null 2>&1; if [ $$? -eq 0 ]; then js-beautify sdr.js/sdr.js >sdr.js/sdr.js.beautiful; mv sdr.js/sdr.js.beautiful sdr.js/sdr.js; fi'

python: libcsdr.so clean-python
	@echo "/* This file is autogenerated. Any changes should be made to csdrpy.i.in */" > csdrpy.i
	@echo 'generating "csdrpy.i"'
	@cat csdrpy.i.in >> csdrpy.i
	@cat fastddc.h >> csdrpy.i
	@cat fft_fftw.h >> csdrpy.i
	@cat fft_rpi.h >> csdrpy.i
	@cat ima_adpcm.h >> csdrpy.i
	@cat libcsdr_gpl.h >> csdrpy.i
	@cat libcsdr.h >> csdrpy.i
	swig -python csdrpy.i
	gcc -c -fPIC csdrpy_wrap.c -I/usr/include/python$(PYVERSION) -I.
	gcc -shared csdrpy_wrap.o -L. -lcsdr -lfftw3f -lm -o _csdrpy.so

clean-python:
	rm -f csdrpy_wrap.c csdrpy_wrap.o _csdrpy.so csdrpy.py* csdrpy.i

python-install: libcsdr.so
	install -m 0755 csdrpy.py _csdrpy.so $(PREFIX)/lib/python$(PYVERSION)/
