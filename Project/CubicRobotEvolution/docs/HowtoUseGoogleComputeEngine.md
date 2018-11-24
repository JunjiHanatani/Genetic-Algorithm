# Google Compute Engine

## Start instance
1. Login [google compute engine][875a6351].
2. Select project and start instance.

## Connect to instance using gcloud compute ssh
1. Download and install [cloud SDK][517812fb], if necessary.
2. Initial setting.
  1. `gcloud init`
  2. ...
3. Connect. <br>
`gcloud compute --project "PROJECT_NAME" ssh --zone "ZONE_NAME" "INSTANCE_NAME"` <br>
Or, a gcloud command line for connecting into the instance can be obtained from "[GCE][875a6351] > VM instances > Instance > Connect > View gcloud command"


## Transfering files between local and instance

##### Local to GCE
`gcloud compute scp [LOCAL_FILE_PATH] [INSTANCE_NAME]:[REMOTE_FILE_PATH]`
##### GCE to local
`gcloud compute scp [INSTANCE_NAME]:[REMOTE_FILE_PATH] [LOCAL_FILE_PATH]
`

`--recurse`: Upload directories recursively. 

## Option - Install applications
`sudo apt-get update`

##### OpenGL
`sudo apt-get install cmake libx11-dev xorg-dev libglu1-mesa-dev freeglut3-dev libglew1.5 libglew1.5-dev libglu1-mesa libglu1-mesa-dev libgl1-mesa-glx libgl1-mesa-dev`

##### OpenMPI
`sudo apt-get install openmpi-doc openmpi-bin libopenmpi-dev`

Execute MPI.
`mpiexec -np 8 -mca btl ^openib "FILE_NAME"`

[517812fb]: https://cloud.google.com/sdk/ "Cloud SDK"
[875a6351]: https://cloud.google.com/compute/ "google compute engine"
