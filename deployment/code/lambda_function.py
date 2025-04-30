import json
import io
import pandas as pd
import matplotlib.pyplot as plt
import boto3
from typing import Dict, Any

def lambda_handler(event: Dict[str, Any], context: Any) -> Dict[str, Any]:
    """
    AWS Lambda handler for RiskMap Genomics PRS calculation.
    
    Args:
        event: Dictionary containing 'bucket' and 'key' for input file
        context: AWS Lambda context object
    
    Returns:
        Dictionary containing success message and output file keys
    """
    # Initialize S3 client
    s3 = boto3.client('s3')
    
    # Get input parameters
    bucket = event['bucket']
    input_key = event['key']
    
    # Download file from S3
    response = s3.get_object(Bucket=bucket, Key=input_key)
    csv_content = response['Body'].read().decode('utf-8')
    
    # Read CSV into DataFrame
    df = pd.read_csv(io.StringIO(csv_content))
    
    # Calculate PRS
    df['PRS'] = df['beta'] * df['effect_allele_freq']
    
    # Save results to CSV in memory
    csv_buffer = io.StringIO()
    df.to_csv(csv_buffer, index=False)
    
    # Upload results CSV to S3
    output_csv_key = 'output/gene_prs_results.csv'
    s3.put_object(
        Bucket=bucket,
        Key=output_csv_key,
        Body=csv_buffer.getvalue()
    )
    
    # Create PRS histogram
    plt.figure(figsize=(10, 6))
    plt.hist(df['PRS'], bins=50, edgecolor='black')
    plt.title('PRS Distribution')
    plt.xlabel('PRS Score')
    plt.ylabel('Frequency')
    
    # Save plot to memory
    img_buffer = io.BytesIO()
    plt.savefig(img_buffer, format='png')
    img_buffer.seek(0)
    
    # Upload plot to S3
    output_img_key = 'output/prs_distribution.png'
    s3.put_object(
        Bucket=bucket,
        Key=output_img_key,
        Body=img_buffer,
        ContentType='image/png'
    )
    
    # Clean up
    plt.close()
    
    return {
        'statusCode': 200,
        'body': json.dumps({
            'message': 'PRS calculation completed successfully',
            'output_files': {
                'results_csv': output_csv_key,
                'distribution_plot': output_img_key
            }
        })
    } 