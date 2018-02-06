variable "org"                    { }
variable "tag"                    { }
variable "name"                   { }
variable "region"                 { }
variable "environment"            { }
variable "allowed_account_ids"    { }
variable "env_vars"               {
  default = "[]"
}
variable "desired_count"          {
  default = "1"
}
variable "cpu"                    {
  default = "128"
}
variable "memory"                 {
  default = "1024"
}
variable "alb_listener_port"      {
  default = "443"
}
variable "container_port"         {
  default = "3838"
}

## Monarch Provider Configurations
provider "aws" {
  version             = "1.2.0"
  region              = "${var.region}"
  allowed_account_ids = [ "${var.allowed_account_ids}" ]
}

provider "logentries" {
  account_key = "${data.terraform_remote_state.stack.le_account_id}"
}

data "aws_ssm_parameter" "dd_api_key" {
  name = "/platform/datadog_api_key"
}

data "aws_ssm_parameter" "dd_app_key" {
  name = "/platform/datadog_app_key"
}

provider "datadog" {
  api_key = "${data.aws_ssm_parameter.dd_api_key.value}"
  app_key = "${data.aws_ssm_parameter.dd_app_key.value}"
}

##

module "triage" {
  source             = "modules/stack/web-application"

  name               = "${var.name}"
  image              = "${data.terraform_remote_state.stack.monarch_repo_short}/${var.org}/${var.name}"
  version            = "${var.tag}"
  port               = "${var.alb_listener_port}"
  container_port     = "${var.container_port}"
  desired_count      = "${var.desired_count}"
  cpu                = "${var.cpu}"
  memory             = "${var.memory}"
  command            = "tini /usr/local/bin/shiny_server.sh"
  environment        = "${data.terraform_remote_state.stack.environment}"
  cluster            = "${data.terraform_remote_state.stack.cluster}"
  iam_role           = "${data.terraform_remote_state.stack.iam_role}"
  security_groups    = "${data.terraform_remote_state.stack.nih_external_elb}"
  subnet_ids         = "${data.terraform_remote_state.stack.external_subnets}"
  log_bucket         = "${data.terraform_remote_state.stack.log_bucket_id}"
  internal_zone_id   = "${data.terraform_remote_state.stack.internal_zone_id}"
  external_zone_id   = "${data.terraform_remote_state.stack.external_zone_id}"
  ssl_certificate_id = "${data.terraform_remote_state.stack.default_ssl_cert_id}"
  vpc_id             = "${data.terraform_remote_state.stack.vpc_id}"
  env_vars           = "${var.env_vars}"
}

module "remote_state" {
  source = "modules/stack/remote-state"

  name              = "${var.name}"
  region            = "${var.region}"
  environment       = "${var.environment}"
}

data "terraform_remote_state" "stack" {
  backend = "s3"

  config {
    bucket  = "${var.environment}-${var.environment}-niaid-terraform-remote-state"
    key     = "${var.environment}-${var.environment}/terraform.tfstate"
    region  = "${var.region}"
  }
}


/**
 * Outputs.
 */

// The listener_id of the ALB created by/for the web-application

output "alb_listener_id" {
  value = "${module.triage.listener_id}"
}

