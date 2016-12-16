#!/bin/bash
ffmpeg -i elec_img/elec_quad_osc_%03d.png -c:v libx264 -crf 24 -pix_fmt yuv420p -r 30 -framerate 30 -t 00:00:08 movies/elec_test.mp4
