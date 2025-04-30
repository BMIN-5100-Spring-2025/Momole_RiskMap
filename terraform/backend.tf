terraform {
  backend "s3" {
    bucket         = "bmin5100-terraform-state"
    key            = "johanna.momole@pennmedicine.upenn.edu-riskmap-genomics/terraform.tfstate"
    region         = "us-east-1"
    encrypt        = true
  }
} 