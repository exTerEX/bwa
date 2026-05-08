.PHONY: all clean test

all:
	cmake -S . -B build
	cmake --build build -j

clean:
	rm -rf build

test: all
	ctest --test-dir build --output-on-failure
