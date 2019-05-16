all:
	@cmake -DCMAKE_BUILD_TYPE=Release -DBENCHMARK_ENABLE_LTO=true -Bbuild -H. && cd build && make

clean:
	trash build bin
