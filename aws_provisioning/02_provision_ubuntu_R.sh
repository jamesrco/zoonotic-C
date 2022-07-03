# 02_provision_ubuntu_R.sh
# provisioning script for an AWS EC2 instance running Ubuntu
# this script sets things up with git, R, and some common R packages for spatial analysis
# created 30 Jun 2022 by Jamie Collins, jcollins@edf.org; tested under Ubuntu 22.04

# *** note: requires at least 16 GB storage just for the install; more if actually trying
# to store anything

cd ~/

export DEBIAN_FRONTEND=noninteractive # in an attempt to get around the interactive prompts that ask us to restart services
# note that we may need to specify as well each time we call something as root (sudo) since this export applies only to current user 

# first, update apt

sudo DEBIAN_FRONTEND=noninteractive apt update

# install git (just in case), AWS CLI
sudo DEBIAN_FRONTEND=noninteractive apt install -y git
echo "Y" | sudo DEBIAN_FRONTEND=noninteractive apt install -y awscli

# install unzip
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y unzip

# install, configure rclone (allows command-line access to Dropbox, Google Drive, etc.)
curl https://rclone.org/install.sh | sudo bash 
echo "q" | rclone config --config=".rclone.conf"

# install R; from https://www.how2shout.com/linux/how-to-install-r-base-ubuntu-22-04-lts-jammy/

sudo DEBIAN_FRONTEND=noninteractive apt install -y --no-install-recommends software-properties-common dirmngr

wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc

sudo DEBIAN_FRONTEND=noninteractive add-apt-repository -y "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

sudo DEBIAN_FRONTEND=noninteractive apt update

echo "Y" | sudo DEBIAN_FRONTEND=noninteractive apt install -y --no-install-recommends r-base

# install OpenMP (so data.table can take advantage of some parallelization)

echo "Y" | sudo DEBIAN_FRONTEND=noninteractive apt install -y libomp-dev

# install necessary packages for working with some spatial data, including data.table
# some help from https://geocompr.github.io/post/2020/installing-r-spatial-ubuntu/

echo "Y" | sudo DEBIAN_FRONTEND=noninteractive apt install -y libudunits2-dev libgdal-dev libgeos-dev libproj-dev libfontconfig1-dev

echo "Y" | sudo DEBIAN_FRONTEND=noninteractive apt install -y r-base-dev r-cran-sf r-cran-raster r-cran-rjava r-cran-data.table