# 01_clone_privaterepo.sh
# script with necessary code to clone a private git repo to an AWS EC2 instance running Ubuntu
# created 30 Jun 2022 by Jamie Collins, jcollins@edf.org; tested under Ubuntu 22.04

# first, update apt

sudo apt update

# install git (just in case), AWS CLI
sudo DEBIAN_FRONTEND=noninteractive apt install -y git
echo "Y" | sudo DEBIAN_FRONTEND=noninteractive apt install -y awscli

# clone the git repo

# # simplest case, if repo is public
#
# git clone https://github.com/jamesrco/zoonotic-c zoonotic-c

# but ... it's more complex if the repo you're looking to clone is private
# easiest way I found that lends itself to automation: use AWS Secrets Manager to store the SSH key you will use to access your private git repo(s), then allow your EC2 instance to access that key securely
# some instructions here: https://medium.com/@nilouferbustani/securing-ssh-private-keys-using-aws-secrets-manager-6d93537c1037
# *** note that the key should be stored as plain text and the JSON brackets {} removed if necessary, since the interface automatically prepends them

# retrieve private key previously stored in Secrets Manager; change --secret-id and --region as needed

aws secretsmanager get-secret-value --secret-id SSH_key_id_ed25519 --query 'SecretString' --region us-west-2 --output text > ~/.ssh/id_ed25519

sudo chmod 600 ~/.ssh/id_ed25519 # necessary since permissions weren't restrictive enough

cd ~/ # move back to home

ssh -o "StrictHostKeyChecking no" github.com

echo "yes" | git clone git@github.com:jamesrco/zoonotic-c.git zoonotic-c # finally, we can clone the private repo

# # start R (if desired)
#
# sudo -i R

