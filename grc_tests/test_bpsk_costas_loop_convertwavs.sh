#!/bin/bash
sox -r 48k -t f32 -c 2 /s/costas_nco -t wav -e floating-point /s/costas_nco.wav 
sox -r 48k -t f32 -c 1 /s/costas_error -t wav -e floating-point /s/costas_error.wav 
sox -r 48k -t f32 -c 1 /s/costas_dphase -t wav -e floating-point /s/costas_dphase.wav 
sox -r 48k -t f32 -c 2 /s/costas_input -t wav -e floating-point /s/costas_input.wav 
sox -r 48k -t f32 -c 2 /s/costas_output -t wav -e floating-point /s/costas_output.wav 



