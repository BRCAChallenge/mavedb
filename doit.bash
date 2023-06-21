#!/usr/bin/env bash
#/Users/melissacline/.local/share/virtualenvs/cool-seq-tool-EK_XvTof
#. /Users/melissacline/.local/share/virtualenvs/cool-seq-tool-EK_XvTof/bin/activate
#. /Users/melissacline/.local/share/virtualenvs/cool-seq-tool-EK_XvTof/bin/activate

# Go to the directory with DynamoDBLocal_lib running and enter this:
# java -Djava.library.path=./DynamoDBLocal_lib -jar DynamoDBLocal.jar -sharedDb

pushd ~/src/python/cool-seq-tool; pipenv shell; popd
export UTA_VERSION=uta_20210129.pgd.gz
export UTA_PASSWORD=uta
export SEQREPO_ROOT_DIR=/usr/local/share/seqrepo/latest
export TRANSCRIPT_MAPPINGS_PATH=${PWD}/cool_seq_tool/data/transcript_mapping.tsv
export AWS_ACCESS_KEY_ID=fakeMyKeyId
export AWS_SECRET_ACCESS_KEY=fakeSecretAccessKey

 
