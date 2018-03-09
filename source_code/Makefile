prog_avx_final:

%: %.cpp
	clang++ --std=gnu++0x -W -Wall -Wno-sign-compare -Ofast -s -pipe -mmmx -msse -msse2 -msse3 -mavx -mcmodel=medium -o $@ $< -pthread

clean:
	rm -f prog_avx_final
