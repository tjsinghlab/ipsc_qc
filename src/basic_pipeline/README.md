## USAGE

cd pipeline

# Build Docker image
docker build -t my_pipeline_image .

# Run all modules
docker run --rm \
  -v $PWD/inputs:/inputs \
  -v $PWD/ref:/ref \
  -v $PWD:/pipeline \
  my_pipeline_image

# Run only cancer and pacnet modules
docker run --rm \
  -v $PWD/inputs:/inputs \
  -v $PWD/ref:/ref \
  -v $PWD:/pipeline \
  my_pipeline_image --no-myco --no-ekaryo --no-outliers


VCF_FILES=(/inputs/*.vcf)
REF_DIR="/ref"
# Results will be in results/ under your repo folder.
