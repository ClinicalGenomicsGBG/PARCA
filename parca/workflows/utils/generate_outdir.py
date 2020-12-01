import yaml
from datetime import datetime
import secrets

class GenerateOutdir:
    @staticmethod
	def read_yaml(filename):
		with open(filename, 'r') as yamlFile:
			yamlDict = yaml.safe_load(yamlFile)

		return yamlDict

	@staticmethod
	def get_date_and_randomizer():
		now = datetime.now()
		dt_string = now.strftime("%Y-%m-%d-%H-%M")
		randomized_value = secrets.token_hex(10)

		sub_outdir=dt_string+randomized_value

		return(sub_outdir)
