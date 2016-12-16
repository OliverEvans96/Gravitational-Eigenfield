#!/bin/bash
ffmpeg -loop 1 -i elec_img/elec_quad_osc_%03d.png -c:v libx264 -crf 24 -pix_fmt yuv420p -t 00:40:00 movies/elec_test.mp4
