variable "org" {}
variable "tag" {}
variable "name" {}
variable "region" {}
variable "environment" {}
variable "allowed_account_ids" {}

variable "docker_labels" {
  default = "null"
}

variable "desired_count" {
  default = "1"
}

variable "cpu" {
  default = "128"
}

variable "memory" {
  default = "1024"
}

variable "alb_listener_port" {
  default = "443"
}

variable "container_port" {
  default = "3838"
}

module "myapp" {
  source        = "modules/stack/web-application"
  docker_labels = "${var.docker_labels}"
}

module "triage" {
  source        = "modules/stack/web-application"
  docker_labels = "${var.docker_labels}"

  name                         = "${var.name}"
  image                        = "${data.terraform_remote_state.stack.monarch_repo_short}/${var.org}/${var.name}"
  image_version                = "${var.tag}"
  port                         = "${var.alb_listener_port}"
  container_port               = "${var.container_port}"
  desired_count                = "${var.desired_count}"
  cpu                          = "${var.cpu}"
  memory                       = "${var.memory}"
  command                      = "tini /usr/local/bin/shiny_server.sh"
  healtcheck_timeout           = "25"
  healtcheck_healthy_threshold = "2"
  environment                  = "${data.terraform_remote_state.stack.environment}"
  cluster                      = "${data.terraform_remote_state.stack.cluster}"
  iam_role                     = "${data.terraform_remote_state.stack.iam_role}"
  security_groups              = "${data.terraform_remote_state.stack.nih_external_elb}"
  subnet_ids                   = "${data.terraform_remote_state.stack.external_subnets}"
  log_bucket                   = "${data.terraform_remote_state.stack.log_bucket_id}"
  internal_zone_id             = "${data.terraform_remote_state.stack.internal_zone_id}"
  external_zone_id             = "${data.terraform_remote_state.stack.external_zone_id}"
  ssl_certificate_id           = "${data.terraform_remote_state.stack.default_ssl_cert_id}"
  vpc_id                       = "${data.terraform_remote_state.stack.vpc_id}"
}

module "remote_state" {
  source = "modules/stack/remote-state"

  name        = "${var.name}"
  region      = "${var.region}"
  environment = "${var.environment}"
}

data "terraform_remote_state" "stack" {
  backend = "s3"

  config {
    bucket = "${var.environment}-${var.environment}-niaid-terraform-remote-state"
    key    = "${var.environment}-${var.environment}/terraform.tfstate"
    region = "${var.region}"
  }
}

/**
 * Outputs.
 */

// The listener_id of the ALB created by/for the web-application

output "alb_listener_id" {
  value = "${module.triage.listener_id}"
}
