"""This script is to draft out functions, classes and general structure of the main CAZy webscraper entry point."""

import sys

import json
import mechanicalsoup


class Site:
	"""A single website, parent to multiple Pages."""
	pages = set()

	def __init__(self, base_url):
		self.base_url = base_url

	def add_page(self, page):
		"""Add (web)page to this (web)Site.

		:param page: Page instancce to add site
		"""
		if page.url not in self.page_urls:
			# do not add duplicate links
			self.pages.add(page)

	def get_class(self, cazy_class):
		"""Return ClassPage for passed CAZy class from this Site.

		ClassPage is the summary page of a CAZy class.

		:param cazy_class: str, name of CAZy class.
		"""
		for page in self.class_pages:
			if page.cazy_class == cazy_class:
				# once the correct Class Page is found
				return page
	def __repr__(self):
		"""Create represent instance."""
		return f"<Site: {id(self)}, base_url: {self.base_url}, page_count: {self.page_count}"

	@property
	def class_pages(self):
		"""Return ClassPages for this site."""
		return [_ for _ in self.pages if isinstance(_, ClassPage)]

	@property
	def page_urls(self):
		"""Return list of all URLs for each page in this Site."""
		return [_.url for _ in self.pages]

	@property
	def page_count(self):
		"""Return the total number of (web)pages in this (web)site."""
		return len(self.pages)

class Page:
	"""Describes the relevant features of a (web)page belonging to the parent (web)Site."""

	links = set()
	url: str

	def __init__(self, url):
		self.url = url

	def add_link(self, url):
		"""Add the passed url to a collection of links retrieved from the current page.

	:param url: str, URL of the link to be added to the collection.
	"""
	self.links.add(url)

	def __repr__(self):
		"""Create represent instance."""
		return f"<Page: {id(self)}, URl={self.url}>"


class ClassPage(Page):
	"""CAZy website page describing and summarising a CAZy class."""

	cazy_class: str

	def __init__(self, url, cazy_class):
		self.cazy_class = cazy_class
		Page.__init(self, url)

	def __repr__(self):
		"""Create represent instance."""
 		return f"<ClassPage: {id(self)}, CLASS={self.cazy_class}, URL={self.url}>"



def main():
	"""Coordinate scraping of CAZy website."""
	# page to start browser: the CAZy homeage
	base_url = "http://www.cazy.org"

	# create browser object
	browser = mechanicalsoupl.Browser()

	# create site class object
	site = Site(base_url)

	# populate site with pages from CAZy homepage
	site = get_links_from_homepage(base_url, browser, site)


def get_links_from_homepage(base_url, browser, site):
	"""Populate site with links from homepage.

	:param base_url: str, url of homepage
	:param browser: bs4 browser object
	:param site: site class object

	Return site object populated with pages from homepage.
	"""

	# scrape CAZy home page
	home_page = retrieve_page(site.base_url, browser, retries, logger)

	# check connection was successul, if not terminate
	if home_page is None:
		logger.error(
			(
				"Failed to connect to CAZy home page.\n"
				"Check network connection.\n"
				"Terminating."
			)
		)
		sys.exit(1)

	# retrieve all links from the home page
	all_links = home_page.soup.select("a")

	# create empty list to store all links to other pages in home page
	pages = []
	for link in all_links:
		try:
			url = link["href"]
			classes = ["Glycoside-Hydrolases.html",
			       "GlycosylTransferases.html",
			       "Polysaccharide-Lyases.html",
			       "Carbohydrate-Esterases.html",
			       "Auxiliary-Activities.html",
			       "Carbohydrate-Binding-Modules.html"]
			if url in classes:
				# if the url is a CAZy class summary page
				# create Class Page class object
				pages.append(ClassPage(url, link.text))
			else:
				# if the url is not for a CAZy class summary page
				# add the url to the list of all pages as a general Page class object
				pages.append(Page(url))
		except KeyError:
			pass

	# populate site object with (web)page objects
	for page in pages:
		site.add_page(page)

	return site


def retrieve_page(url, browser, retries, logger):
	"""Try connection to url, retry if fails, retrieve page.

	:param url: str, url for webpage
	:param browser: bs4 browser object
	:param args: parser obejct

	Return response object.
	"""
	# create reponse object
	page = browser.get(url)

	# check connection was successfully made
	# <Response [200]> marks success
	if str(page) != "<Response [200]>":
		tries = 0
		while (str(page) != "<Response [200]>") and (tries < retries):
			page = browser.get(url)
			tries += 1

	if str(page) != "<Response [200]>":
		logger.error(
			(
				"Network error encounted too many times when trying to connect to,\n"
				f"{url}"
			)
		)
		return None
	
	return page

def parse_config_file(config, logger):
	## create function to parse config and specify classes and families
