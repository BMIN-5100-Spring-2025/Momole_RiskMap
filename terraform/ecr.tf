resource "aws_ecr_repository" "riskmap_genomics" {
  name = "johanna_riskmap_genomics"
  
  image_scanning_configuration {
    scan_on_push = true
  }

  tags = {
    Name        = "johanna_riskmap_genomics"
    Environment = var.environment
    Owner       = "johanna"
  }
}

resource "aws_ecr_lifecycle_policy" "riskmap_genomics" {
  repository = aws_ecr_repository.riskmap_genomics.name

  policy = jsonencode({
    rules = [{
      rulePriority = 1
      description  = "Keep last 30 images"
      selection = {
        tagStatus   = "any"
        countType   = "imageCountMoreThan"
        countNumber = 30
      }
      action = {
        type = "expire"
      }
    }]
  })
} 