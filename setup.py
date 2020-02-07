
import site
import subprocess

# build with cargo
ret = subprocess.call(["cargo", "build", "--release"])

if ret == 0:
    # Insert into site package
    location = site.getsitepackages()[0]
    subprocess.call(["sudo", "mv", "target/release/libdenoise_engine.so", "{}/denoise_engine.so".format(location)])

    
