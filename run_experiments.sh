EXPERIMENTS=experiments/
IMAGE_DIMENSIONS=(50 75 100)

function scale() {
    dim=$1
    mkdir -p cropped/{bdr,bdc}
    python crop_images.py data/bdc data/coords.3368.txt cropped/bdc ${dim}x${dim}
    python crop_images.py data/bdr data/coords.3368.txt cropped/bdr ${dim}x${dim}
}

for dim in ${IMAGE_DIMENSIONS[@]}; do
    echo "Running experiments for images of size $dim"
    resultsd=$EXPERIMENTS/$dim
    mkdir -p $resultsd

    rm -rf cropped
    scale $dim
    r --save < face-recognition.r | tail -n 10 > $resultsd/results
    mv .RData $resultsd/RData
done

rm -rf cropped
