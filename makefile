all: release

release:
	@cmake -DCMAKE_BUILD_TYPE=Release -DBENCHMARK_ENABLE_LTO=true -Bbuild -H. -DCMAKE_EXPORT_COMPILE_COMMANDS=YES && cd build && make

mpi:
	@cmake -DMPI_BUILD=ON -DCMAKE_BUILD_TYPE=Release -DBENCHMARK_ENABLE_LTO=true -Bbuild -H. -DCMAKE_EXPORT_COMPILE_COMMANDS=YES && cd build && make

debug:
	@cmake -DCMAKE_BUILD_TYPE=Debug -DBENCHMARK_ENABLE_LTO=true -Bbuild -H. -DCMAKE_EXPORT_COMPILE_COMMANDS=YES && cd build && make

static-gsl:
	@cmake -DCMAKE_BUILD_TYPE=Release -DSTATIC_GSL=On -DBENCHMARK_ENABLE_LTO=true -Bbuild -H. -DCMAKE_EXPORT_COMPILE_COMMANDS=YES && cd build && make

static:
	@cmake -DCMAKE_BUILD_TYPE=Release -DSTATIC_BUILD=On  -Bbuild -H. -DCMAKE_EXPORT_COMPILE_COMMANDS=YES && cd build && make


clean:
	rm -rf build bin
