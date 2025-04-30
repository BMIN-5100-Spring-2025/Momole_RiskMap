# RiskMap Genomics AWS Lambda Function

This Lambda function processes GWAS data to calculate Polygenic Risk Scores (PRS) and generate visualizations.

## Functionality

- Accepts S3 bucket and key for input GWAS CSV file
- Calculates PRS using beta and effect allele frequency
- Generates a PRS distribution histogram
- Saves results back to S3

## Deployment Instructions

1. Create a deployment package:
   ```bash
   # Create a virtual environment
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate

   # Install dependencies
   pip install -r requirements.txt

   # Create deployment package
   pip install -r requirements.txt --target deployment/package
   cp lambda_function.py deployment/package/
   cd deployment/package
   zip -r ../riskmap_genomics.zip .
   ```

2. Create an S3 bucket for input/output files:
   - Create a bucket named `riskmap-bucket` (or your preferred name)
   - Create `input/` and `output/` folders in the bucket

3. Create the Lambda function:
   - Runtime: Python 3.10
   - Memory: 1024MB
   - Timeout: 1 minute
   - Upload the `riskmap_genomics.zip` deployment package
   - Set the handler to `lambda_function.lambda_handler`

4. Configure IAM role:
   - Create an IAM role with S3 read/write permissions
   - Attach the role to the Lambda function

## Testing

Test the function with the following JSON payload:
```json
{
    "bucket": "riskmap-bucket",
    "key": "input/sample.csv"
}
```

## Input File Format

The input CSV should contain at least these columns:
- `beta`: Effect size
- `effect_allele_freq`: Effect allele frequency

## Output

The function generates:
1. `output/gene_prs_results.csv`: CSV with PRS calculations
2. `output/prs_distribution.png`: Histogram of PRS distribution
