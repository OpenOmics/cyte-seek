## Steps for Building Docker Images

Directly below are instructions for building a base image for the `actinn cite` pipeline using the provided 
Dockerfile:

```bash
# See listing of images on computer
docker image ls

# Build from Dockerfile
docker build --no-cache -f Dockerfile --tag=actinn:v0.1.0 .

# Testing, take a peek inside
docker run -ti actinn:v0.1.0 /bin/bash

# Updating Tag  before pushing to DockerHub
docker tag actinn:v0.1.0 skchronicles/actinn:v0.1.0
docker tag actinn:v0.1.0 skchronicles/actinn         # latest

# Check out new tag(s)
docker image ls

# Push new tagged image to DockerHub
docker push skchronicles/actinn:v0.1.0
docker push skchronicles/actinn:latest
```

### Other Recommended Steps

Scan your image for known vulnerabilities:

```bash
docker scan actinn:v0.1.0
```

> **Please Note**: Any references to `skchronicles` should be replaced your username if you would also like to push the image to a non-org account.
