all: lcm.c
	gcc -O3 -mavx -o topkwy-dfs -DLCMFREQ_COUNT_OUTPUT lcm.c -lm
profile: lcm.c
	gcc -O3 -mavx -o topkwy-dfs -DLCMFREQ_COUNT_OUTPUT -DPROFILE_MINING lcm.c -lm
