all:
	@cmake -DCMAKE_BUILD_TYPE=Release -DBENCHMARK_ENABLE_LTO=true -Bbuild -H. && cd build && make

debug:
	@cmake -DCMAKE_BUILD_TYPE=Debug -DBENCHMARK_ENABLE_LTO=true -Bbuild -H. && cd build && make

static-gsl:
	@cmake -DCMAKE_BUILD_TYPE=Release -DSTATIC_GSL=On -DBENCHMARK_ENABLE_LTO=true -Bbuild -H. && cd build && make

clean:
	trash build bin
