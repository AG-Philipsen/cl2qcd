/** @file
 * Definition of utilities to normalize input file contents
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#ifndef _CONFIG_FILE_NORMALIZER_
#define _CONFIG_FILE_NORMALIZER_

#include <map>
#include <istream>

namespace meta {
	/**
	 * This class provides utilities to normalize the input file.
	 *
	 * This is mostly a workaround for boost::parse_config_file not allowing additional parsers.
	 * Therefore the normalization only applies to the option names, not the values.
	 *
	 * To use first define a set of aliases on an instance of ConfigFileNormalizer.
	 *
	 * \code
	 * ConfigFileNormalizer normalizer;
	 * normalizer.add_alias("foobar", "foo");
	 * normalizer.add_alias("bar", "foo");
	 * \endcode
	 *
	 * Then invoke the ConfigFileNormalizer instance on some file to get the normalized input stream.
	 * The normalizer will replace all aliases by their mapped values.
	 *
	 * \code
	 * std::istream normalized_stream = normalizer("some_file");
	 * \endcode
	 */
	class ConfigFileNormalizer {
		public:
			/**
			 * Create an instance that does not yet contain any alias information.
			 */
			ConfigFileNormalizer() : aliases() { };

			/**
			 * Add an alias to the ConfigFileNormalizer instance.
			 *
			 * When the ConfigFileNormalizer is invoked it will replace all instances of the
			 * aliased option by the mapped one.
			 *
			 * \param alias The option to replace
			 * \param mapped The value to replace the option by
			 */
			void add_alias(std::string alias, std::string mapped);

			/**
			 * Get the normalized contents of the given file.
			 *
			 * \param file The name of the file to retrieve.
			 * \return A string containing the normalized config file
			 *
			 * \todo add case normalization
			 */
			std::string operator() (std::string file) const;

			/**
			 * Get the normalized contents of the given file.
			 *
			 * \param input Stream of the contents to normalize
			 * \return A string containing the normalized config file
			 *
			 * \todo add case normalization
			 */
			std::string operator() (std::istream& input) const;

		private:
			/**
			 * Internal storage of the alias names with their mapped values
			 */
			std::map<std::string, std::string> aliases;
	};
}

#endif /* _CONFIG_FILE_NORMALIZER_ */
