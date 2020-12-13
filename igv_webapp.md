## IGV-webapp Installation

### Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh

### create a new conda environment (`igv`) and install required packages

    conda create -n igv
    conda activate igv
    conda install nodejs
    conda list

### Clone the igv-webapp repository to local:

    cd /home/springer/zhoux379/git
    git clone git@github.com:igvteam/igv-webapp.git

(or if you cloned the repo a while ago) update the local repo (if there are new changes):

    cd ./igv-webapp
    git pull

### Build and install

    cd ./igv-webapp
    npm install
    npm run build

### Create an empty S3 bucket served on msi.umn.edu


    # make a new bucket
    s3cmd mb s3://zhoup-igv-test

    # inside the `igv-webapp` directory
    s3cmd sync --delete-removed --acl-public --follow-symlinks --exclude '.git/* .github/*' --no-mime-magic --guess-mime-type dist/ s3://zhoup-igv-test

    # check if files have been successfully uploaded
    s3cmd ls s3://zhoup-igv-test

    # double check everything is public
    s3cmd setacl --acl-public -r s3://zhoup-igv-test

Now the app should be accessible online:
	https://s3.msi.umn.edu/zhoup-igv-test/index.html

### Customized genome configuration
# list your existing buckets
s3cmd ls
- Check [igv-webapp github](https://github.com/igvteam/igv-webapp)
- Here is [a sample configuration for maize](https://github.com/orionzhou/igv-webapp/blob/master/igvwebConfig.js)
