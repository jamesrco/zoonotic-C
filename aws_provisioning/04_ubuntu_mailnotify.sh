# 04_ubuntu_mailnotify.sh

# provisioning script for an AWS EC2 instance running Ubuntu
# this script installs necessary tools & files to allow EC2 instance to send me
# notification emails

# created 1 Jul 2022 by Jamie Collins, jcollins@edf.org; tested under Ubuntu 22.04
# with help from https://www.r-bloggers.com/2018/09/tired-of-waiting-for-your-r-scripts-to-finish-let-aws-do-the-work-get-notified-by-e-mail/
# and https://gist.github.com/titpetric/114eb27f6e453e3e8849d65ca1a3d360
# with modification to use an app password since Google has officially ended support
# for the "less secure email access" option

# *** in any case, use a throwaway email address and a "one time" or app password
# since password has to be send via plain text

# install ssmtp

sudo DEBIAN_FRONTEND=noninteractive apt install -y ssmtp

# create ssmtp.conf file

# set values

SSMTP_ROOT=${SSMTP_ROOT:-jamesrco.aws@gmail.com}
SSMTP_USER=${SSMTP_USER:-jamesrco.aws@gmail.com}
SSMTP_PASS=${SSMTP_PASS:-***REMOVED***}
SSMTP_SERVER=${SSMTP_SERVER:-smtp.gmail.com}
SSMTP_PORT=${SSMTP_PORT:-587}
SSMTP_HOSTNAME=${echo \$(hostname -f)}

# concatenate; write file to the right place

sudo cat << EOF > /etc/ssmtp/ssmtp.conf
root=$SSMTP_ROOT
mailhub=$SSMTP_SERVER:$SSMTP_PORT
AuthUser=$SSMTP_USER
AuthPass=$SSMTP_PASS
UseSTARTTLS=YES
UseTLS=YES
hostname=$SSMTP_HOSTNAME
FromLineOverride=YES
EOF

# sudo cp ~/zoonotic-c/aws_provisioning/ssmtp/ssmtp.conf /etc/ssmtp/ssmtp.conf
