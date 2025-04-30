resource "aws_s3_bucket" "momole_riskmap_data" {
  bucket = "momole-riskmap-data"

  tags = {
    Owner = element(split("/", data.aws_caller_identity.current.arn), 1)
  }
}

resource "aws_s3_bucket_versioning" "momole_riskmap_data" {
  bucket = aws_s3_bucket.momole_riskmap_data.id
  versioning_configuration {
    status = "Enabled"
  }
}

resource "aws_s3_bucket_server_side_encryption_configuration" "momole_riskmap_data" {
  bucket = aws_s3_bucket.momole_riskmap_data.id

  rule {
    apply_server_side_encryption_by_default {
      sse_algorithm = "AES256"
    }
  }
}

# Import remote state for shared infrastructure
data "terraform_remote_state" "infrastructure" {
  backend = "s3"
  config = {
    bucket = "bmin5100-terraform-state"
    key    = "infrastructure/terraform.tfstate"
    region = "us-east-1"
  }
}

# CloudWatch Log Group
resource "aws_cloudwatch_log_group" "riskmap_genomics" {
  name              = "/ecs/riskmap-genomics"
  retention_in_days = 30
}

# ECS Task Definition
resource "aws_ecs_task_definition" "riskmap_genomics" {
  family                   = "johanna_riskmap_genomics"
  network_mode             = "awsvpc"
  requires_compatibilities = ["FARGATE"]
  cpu                      = var.container_cpu
  memory                   = var.container_memory
  execution_role_arn       = aws_iam_role.ecs_task_execution_role.arn
  task_role_arn            = aws_iam_role.ecs_task_role.arn

  ephemeral_storage {
    size_in_gib = 100
  }

  container_definitions = jsonencode([
    {
      name  = "riskmap-genomics"
      image = "${aws_ecr_repository.riskmap_genomics.repository_url}:latest"
      environment = [
        {
          name  = "S3_BUCKET"
          value = "arn:aws:s3:::${var.s3_bucket_name}"
        }
      ]
      logConfiguration = {
        logDriver = "awslogs"
        options = {
          awslogs-group         = aws_cloudwatch_log_group.riskmap_genomics.name
          awslogs-region        = var.aws_region
          awslogs-stream-prefix = "ecs"
        }
      }
    }
  ])
}

# ECS Service
resource "aws_ecs_service" "riskmap_genomics" {
  name            = "johanna_riskmap_genomics"
  cluster         = data.terraform_remote_state.infrastructure.outputs.ecs_cluster_id
  task_definition = aws_ecs_task_definition.riskmap_genomics.arn
  desired_count   = 1
  launch_type     = "FARGATE"

  network_configuration {
    subnets          = [data.terraform_remote_state.infrastructure.outputs.private_subnet_id]
    security_groups  = [data.terraform_remote_state.infrastructure.outputs.ecs_security_group_id]
    assign_public_ip = false
  }

  depends_on = [
    aws_iam_role_policy_attachment.ecs_task_execution_role_policy,
    aws_iam_role_policy.s3_access_policy
  ]
} 