variable "aws_region" {
  description = "AWS region"
  type        = string
  default     = "us-east-1"
}

variable "environment" {
  description = "Environment name"
  type        = string
  default     = "production"
}

variable "app_name" {
  description = "Application name"
  type        = string
  default     = "riskmap-genomics"
}

variable "aws_account_id" {
  description = "AWS account ID"
  type        = string
}

variable "s3_bucket_name" {
  description = "S3 bucket name for data storage"
  type        = string
  default     = "riskmap-genomics-johanna"
}

variable "container_cpu" {
  description = "Container CPU units"
  type        = number
  default     = 512  # 0.5 vCPU
}

variable "container_memory" {
  description = "Container memory in MB"
  type        = number
  default     = 1024  # 1GB
}

variable "container_port" {
  description = "Container port"
  type        = number
  default     = 80
} 