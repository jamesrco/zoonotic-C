# 03_ubuntu_swapmake.sh

# provisioning script for an AWS EC2 instance running Ubuntu
# this script creates a swapfile to dedicate disk space as memory
# useful to avert killing of R sessions due to out of memory errors

# created 1 Jul 2022 by Jamie Collins, jcollins@edf.org; tested under Ubuntu 22.04
# with help from https://www.digitalocean.com/community/tutorials/how-to-add-swap-space-on-ubuntu-20-04

# assumes you have at least as much storage on your image as dedicated to a swap

sudo fallocate -l 35G /swapfile

sudo chmod 600 /swapfile

sudo mkswap /swapfile

sudo swapon /swapfile

# to make sure it's enabled and working:

sudo swapon --show

free -h