# 02_provision_ubuntu_R.sh
# provisioning script for an AWS EC2 instance running Ubuntu
# this script sets things up with git, R, and some common R packages for spatial analysis
# created 30 Jun 2022 by Jamie Collins, jcollins@edf.org; tested under Ubuntu 22.04

# first, update apt

sudo apt update

# install git (just in case), AWS CLI
sudo apt install git
echo "Y" | sudo apt install awscli

# install unzip
sudo apt-get install unzip

# install, configure rclone (allows command-line access to Dropbox, Google Drive, etc.)
curl https://rclone.org/install.sh | sudo bash 
echo "q" | rclone config --config=".rclone.conf"

# install R; from https://www.digitalocean.com/community/tutorials/how-to-install-r-on-ubuntu-22-04

export DEBIAN_FRONTEND=noninteractive # so we can get around the interactive prompts that ask us to restart services

wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo gpg --dearmor -o /usr/share/keyrings/r-project.gpg

echo "deb [signed-by=/usr/share/keyrings/r-project.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" | sudo tee -a /etc/apt/sources.list.d/r-project.list

echo "Y" | sudo apt install --no-install-recommends r-base

# install necessary packages for working with some spatial data, including data.table
# some help from https://geocompr.github.io/post/2020/installing-r-spatial-ubuntu/

echo "Y" | sudo apt install libudunits2-dev libgdal-dev libgeos-dev libproj-dev libfontconfig1-dev

echo "Y" | sudo apt install r-base-dev r-cran-sf r-cran-raster r-cran-rjava r-cran-data.table