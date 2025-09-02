# docker execution
bash docker_execute.sh 24

# stop container
docker container stop 20250807_Fangchao

# view docker log:
sudo vi `docker inspect --format='{{.LogPath}}' 20250807_Fangchao`

# SSH port forwarding
ssh -L 9090:localhost:9090 shiyu@192.168.11.6

# build image
docker build . -t multi_organ_mfa --build-arg HTTP_PROXY="http://127.0.0.1:7897" \
  --build-arg HTTPS_PROXY="http://127.0.0.1:7897"
