
#!/bin/bash
lefse_format_input.py input_lefse.txt input_lefse.in -c 1 -u 2 -o 1000000
lefse_run.py input_lefse.in lefse.res -a 0.1 -w 0.1
