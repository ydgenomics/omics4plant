
```shell
docker pull peevin/spotclean:v1
docker login -u ydgenomics
# 查看镜像是否装了jupyter
# 注意：这会启动一个临时container，执行完命令后会自动退出并删除（没有-d参数）
docker run peevin/spotclean:v1 jupyter --version

# 进入container进行操作
docker run -it peevin/spotclean:v1 /bin/bash
# 在container内执行操作，例如：
# sudo apt update && sudo apt install -y python3-pip r-base r-base-dev && python3 -m pip install --upgrade notebook jupyterlab && R -e "install.packages('IRkernel', repos='https://cloud.r-project.org'); IRkernel::installspec(user = FALSE)"
# jupyter notebook --ip=0.0.0.0 --allow-root
# 安装完成后，退出container（exit或Ctrl+D）

# 保存修改后的image
# 先查看container ID
docker ps -a
# 使用commit命令保存为新image
docker commit 58c345390814 ydgenomics/spotclean:v1
docker run ydgenomics/spotclean:v1 jupyter -v

# 推送镜像到DockerHub
# 先给镜像打标签（格式：username/repository:tag）
# docker tag的意义：为镜像创建别名，指向DockerHub的用户/仓库位置，以便push时知道推送到哪里
# 必须加：push时Docker会根据tag中的username识别推送目标，不加tag就无法推送到DockerHub
docker tag ydgenomics/spotclean:v1 ydgenomics/spotclean:v1
docker push ydgenomics/spotclean:v1
```

image是静态的，container是动态的

```shell
sudo systemctl status docker

sudo nano /etc/docker/daemon.json
{
 "insecure-registries": ["registry.example.com"]
}
sudo systemctl daemon-reload
sudo systemctl restart docker
```

- 一个文件搞定80%的Docker生产问题：daemon.json配置全解析 https://mp.weixin.qq.com/s/lze5nAsMVz8k4vYQym1XJQ