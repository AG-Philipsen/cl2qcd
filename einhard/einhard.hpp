/**
 * @file
 *
 * This is the main include file for Einhard.
 *
 * Copyright 2010 Matthias Bach <marix@marix.org>
 *
 * This file is part of Einhard.
 *
 * Einhard is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Einhard is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Einhard.  If not, see <http://www.gnu.org/licenses/>.
 */

/** \mainpage Einhard Logging Library
 *
 * \section intro_sec Introduction
 *
 * Einhard aims at being a lightweight logging library with the following features:
 *  - Severity filtering
 *  - Automatic output colorization
 *  - Output timestamping
 *  - Low performance overhead
 *
 * To use Einhard all you need to do is create an einhard::Logger object and push messages into its output
 * streams.
 *
 * \code
 * using namespace einhard;
 * Logger logger( INFO );
 * logger.trace() << "Trace message"; // will not be printed
 * logger.error() << "Error message"; // will be printed
 * \endcode
 *
 * To reduce the performance overhad you can specify NDEBUG during compile. This will disable trace
 * and debug messages in a way that should allow the compilers dead code elimination to remove
 * everything that is only pushed into such a stream.
 *
 * You can use all of the normal C++ stream modifiers like std::setw, std::setfill and std::setprecision.
 *
 * \code
 * using namespace std;
 * ...
 * logger.info() << "Iterations: " << setfill(' ') << setw(5) << count;
 * logger.info() << "Load: " << setprecision(1) << load << "%";
 * \endcode
 *
 * \section install_sec Installation
 *
 * Einhard is build using cmake. You can install it using the usual cmake triplet:
 * |code
 * cmake
 * make
 * sudo make install
 * |endcode
 * If you want to build a static lib or install to a custom path you can use the usual cmake
 * utilities to adjust the configuration.
 */

#include <iostream>
#include <iomanip>
#include <ctime>

// This C header is sadly required to check whether writing to a terminal or a file
#include <stdio.h>
#include <unistd.h>

/**
 * This namespace contains all objects required for logging using Einhard.
 */
namespace einhard {
/**
   * Version string of the Enhard library
 */
static char const VERSION[] = "0.2";

/**
   * Specification of the message severity.
   *
   * In output each level automatically includes the higher levels.
   */
enum LogLevel { ALL, /**< Log all message */
                TRACE, /**< The lowes severity for messages describing the program flow */
                DEBUG, /**< Debug messages */
                INFO, /**< Messages of informational nature, expected processing time e.g. */
                WARN, /**< Warning messages */
                ERROR, /**< Non-fatal errors */
                FATAL, /**< Messages that indicate terminal application failure */
                OFF /**< If selected no messages will be output */
              };

/**
   * Retrieve a human readable representation of the given log level value.
   */
inline char const * getLogLevelString( const LogLevel level );

/**
   * A stream modifier that allows to colorize the log output.
   */
template<typename Parent> struct Color {
	bool reset;
	Color() : reset(true) {}
	Color<Parent> operator~() const {
		Color<Parent> tmp(*this);
		tmp.reset = false;
		return tmp;
	}
	char const *ansiCode() const {
		return Parent::ANSI;
	}
	bool resetColor() const {
		return reset;
	}
};
#define _COLOR(name, code) \
  struct _##name { static char const ANSI[]; }; \
  typedef Color<_##name> name

_COLOR(DGray,   "01;30");
_COLOR(Black,   "00;30");
_COLOR(Red,     "01;31");
_COLOR(DRed,    "00;31");
_COLOR(Green,   "01;32");
_COLOR(DGreen,  "00;32");
_COLOR(Yellow,  "01;33");
_COLOR(Orange,  "00;33");
_COLOR(Blue,    "01;34");
_COLOR(DBlue,   "00;34");
_COLOR(Magenta, "01;35");
_COLOR(DMagenta, "00;35");
_COLOR(Cyan,    "01;36");
_COLOR(DCyan,   "00;36");
_COLOR(White,   "01;37");
_COLOR(Gray,    "00;37");
_COLOR(NoColor, "0"    );
#undef _COLOR

#define ANSI_COLOR_WARN  _Orange::ANSI
#define ANSI_COLOR_ERROR _DRed::ANSI
#define ANSI_COLOR_FATAL _DRed::ANSI
#define ANSI_COLOR_INFO  _DGreen::ANSI
#define ANSI_COLOR_DEBUG _DBlue::ANSI
#define ANSI_COLOR_CLEAR _NoColor::ANSI

/**
 * A wrapper for the output stream taking care proper formatting and colorization of the output.
 */
template<LogLevel VERBOSITY> class OutputFormatter {
private:
	// The output stream to print to
	std::ostream * const out;
	// Whether to colorize the output
	bool const colorize;
	mutable bool resetColor;

public:
	OutputFormatter( std::ostream * const out, bool const colorize ) : out( out ), colorize( colorize ),
		resetColor(false) {
		if( out != 0 ) {
			// Figure out current time
			time_t rawtime;
			tm * timeinfo;
			time( &rawtime );
			timeinfo = localtime( &rawtime );

			if( colorize ) {
				// set color according to log level
				switch ( VERBOSITY ) {
					case WARN:
						*out << ANSI_COLOR_WARN;
						break;
					case ERROR:
						*out << ANSI_COLOR_ERROR;
						break;
					case FATAL:
						*out << ANSI_COLOR_FATAL;
						break;
					case INFO:
						*out << ANSI_COLOR_INFO;
						break;
					case DEBUG:
						*out << ANSI_COLOR_DEBUG;
						break;
					default:
						// in other cases we leave the default color
						break;
				}
			}

			// output it
			*out << '[';
			*out << std::setfill('0') << std::setw(2) << timeinfo->tm_hour;
			*out << ':';
			*out << std::setfill('0') << std::setw(2) << timeinfo->tm_min;
			*out << ':';
			*out << std::setfill('0') << std::setw(2) << timeinfo->tm_sec;
			*out << ']';
			// TODO would be good to have this at least .01 seconds
			// for non-console output pure timestamp would probably be better

			// output the log level of the message
			*out << ' ' << getLogLevelString( VERBOSITY ) << ": ";

			if( colorize ) {
				// reset color according to log level
				switch ( VERBOSITY ) {
					case INFO:
					case DEBUG:
					case WARN:
					case ERROR:
					case FATAL:
						*out << ANSI_COLOR_CLEAR;
						break;
					default:
						// in other cases color is still default anyways
						break;
				}
			}
		}
	}

	template<typename T> inline const OutputFormatter<VERBOSITY>& operator<<( const Color<T> &col ) const {
		if (out && colorize) {
			*out << col.ansiCode();
			resetColor = col.resetColor();
		}
		return *this;
	}

	template<typename T> inline const OutputFormatter<VERBOSITY>& operator<<( const T &msg ) const {
		// output the log message
		if( out != 0 ) {
			*out << msg;
			if (resetColor) {
				*out << ANSI_COLOR_CLEAR;
				resetColor = false;
			}
		}

		return *this;
	}

	~OutputFormatter( ) {
		if( out != 0 ) {
			// make sure there is no strange formatting set anymore
			*out << std::resetiosflags(  std::ios_base::floatfield  | std::ios_base::basefield
			                             | std::ios_base::adjustfield | std::ios_base::uppercase
			                             | std::ios_base::showpos     | std::ios_base::showpoint
			                             | std::ios_base::showbase    |  std::ios_base::boolalpha );

			*out << std::endl; // TODO this would probably better be only '\n' as to not flush the buffers
		}
	}
};

/**
   * A Logger object can be used to output messages to stdout.
   *
   * The Logger object is created with a certain verbosity. All messages of a more verbose
   * LogLevel will be filtered out. The way the class is build this can happen at compile
   * time up to the level restriction given by the template parameter.
   *
   * The class can automatically detect non-tty output and will not colorize output in that case.
   */
template < LogLevel MAX = ALL > class Logger {
private:
	LogLevel verbosity;
	bool colorize;

public:
	/**
	 * Create a new Logger object.
	 *
	 * The object will automatically colorize output on ttys and not colorize output
	 * on non ttys.
	 */
	Logger( const LogLevel verbosity = WARN ) : verbosity( verbosity ) {
		// use some, sadly not c++-ways to figure out whether we are writing ot a terminal
		// only colorize when we are writing ot a terminal
		colorize = isatty( fileno( stdout ) );
	};
	/**
	 * Create a new Logger object explicitly selecting whether to colorize the output or not.
	 *
	 * Be aware that if output colorization is selected output will even be colorized if
	 * output is to a non tty.
	 */
	Logger( const LogLevel verbosity, const bool colorize ) : verbosity( verbosity ), colorize( colorize ) { };

	/** Access to the trace message stream. */
	inline const OutputFormatter<TRACE> trace() const;
	/** Access to the debug message stream. */
	inline const OutputFormatter<DEBUG> debug() const;
	/** Access to the info message stream. */
	inline const OutputFormatter<INFO> info() const;
	/** Access to the warning message stream. */
	inline const OutputFormatter<WARN> warn() const;
	/** Access to the error message stream. */
	inline const OutputFormatter<ERROR> error() const;
	/** Access to the fatal message stream. */
	inline const OutputFormatter<FATAL> fatal() const;

	/** Check whether trace messages will be output */
	inline bool beTrace() const;
	/** Check whether debug messages will be output */
	inline bool beDebug() const;
	/** Check whether info messages will be output */
	inline bool beInfo() const;
	/** Check whether warn messages will be output */
	inline bool beWarn() const;
	/** Check whether error messages will be output */
	inline bool beError() const;
	/** Check whether fatal messages will be output */
	inline bool beFatal() const;

	/** Modify the verbosity of the Logger.
	 *
	 * Be aware that the verbosity can not be increased over the level given by the template
	 * parameter
	 */
	inline void setVerbosity( LogLevel verbosity ) {
		this->verbosity = verbosity;
	}
	/** Retrieve the current log level.
	 *
	 * If you want to check whether messages of a certain LogLevel will be output the
	 * methos beTrace(), beDebug(), beInfo(), beWarn(), beError() and beFatal() should be
	 * preffered.
	 */
	inline LogLevel getVerbosity( ) const {
		return this->verbosity;
	}
	/**
	 * Retrieve a human readable representation of the current log level
	 */
	inline char const * getVerbosityString( ) const {
		return getLogLevelString( this->verbosity );
	}
	/**
	 * Select whether the output stream should be colorized.
	 */
	inline void setColorize( bool colorize ) {
		this->colorize = colorize;
	}
	/**
	 * Check whether the output stream is colorized.
	 */
	inline bool getColorize( ) const {
		return this->colorize;
	}
};

/*
 * IMPLEMENTATIONS
 */

inline char const * getLogLevelString( const LogLevel level )
{
	switch( level ) {
		case ALL:
			return "ALL";
		case TRACE:
			return "TRACE";
		case DEBUG:
			return "DEBUG";
		case INFO:
			return "INFO";
		case WARN:
			return "WARNING";
		case ERROR:
			return "ERROR";
		case FATAL:
			return "FATAL";
		case OFF:
			return "OFF";
		default:
			return "--- Something is horribly broken - A value outside the enumation has been given ---";
	}
}

template<LogLevel MAX> const OutputFormatter<TRACE> Logger<MAX>::trace() const
{
	if( beTrace() )
		return OutputFormatter<TRACE>( &std::cout, colorize );
	else
		return OutputFormatter<TRACE>( 0, colorize );
}

template<LogLevel MAX> bool Logger<MAX>::beTrace() const
{
#ifdef NDEBUG
	return false;
#else
	return ( MAX <= TRACE && verbosity <= TRACE );
#endif
}

template<LogLevel MAX> const OutputFormatter<DEBUG> Logger<MAX>::debug() const
{
	if( beDebug() )
		return OutputFormatter<DEBUG>( &std::cout, colorize );
	else
		return OutputFormatter<DEBUG>( 0, colorize );
}

template<LogLevel MAX> bool Logger<MAX>::beDebug() const
{
#ifdef NDEBUG
	return false;
#else
	return ( MAX <= DEBUG && verbosity <= DEBUG );
#endif
}

template<LogLevel MAX> const OutputFormatter<INFO> Logger<MAX>::info() const
{
	if( beInfo() )
		return OutputFormatter<INFO>( &std::cout, colorize );
	else
		return OutputFormatter<INFO>( 0, colorize );
}

template<LogLevel MAX> bool Logger<MAX>::beInfo() const
{
	return ( MAX <= INFO && verbosity <= INFO );
}

template<LogLevel MAX> const OutputFormatter<WARN> Logger<MAX>::warn() const
{
	if( beWarn() )
		return OutputFormatter<WARN>( &std::cout, colorize );
	else
		return OutputFormatter<WARN>( 0, colorize );
}

template<LogLevel MAX> bool Logger<MAX>::beWarn() const
{
	return ( MAX <= WARN && verbosity <= WARN );
}

template<LogLevel MAX> const OutputFormatter<ERROR> Logger<MAX>::error() const
{
	if( beError() )
		return OutputFormatter<ERROR>( &std::cout, colorize );
	else
		return OutputFormatter<ERROR>( 0, colorize );
}

template<LogLevel MAX> bool Logger<MAX>::beError() const
{
	return ( MAX <= ERROR && verbosity <= ERROR );
}

template<LogLevel MAX> const OutputFormatter<FATAL> Logger<MAX>::fatal() const
{
	if( beFatal() )
		return OutputFormatter<FATAL>( &std::cout, colorize );
	else
		return OutputFormatter<FATAL>( 0, colorize );
}

template<LogLevel MAX> bool Logger<MAX>::beFatal() const
{
	return ( MAX <= FATAL && verbosity <= FATAL );
}

}

#undef ANSI_COLOR_WARN
#undef ANSI_COLOR_ERROR
#undef ANSI_COLOR_FATAL
#undef ANSI_COLOR_INFO
#undef ANSI_COLOR_DEBUG
#undef ANSI_COLOR_CLEAR

#define EINHARD_STREAM(arg) \
template<einhard::LogLevel VERBOSITY> \
inline const einhard::OutputFormatter<VERBOSITY>& operator<<(const einhard::OutputFormatter<VERBOSITY> &out, arg)
#define EINHARD_STREAM_T1(T1, arg) \
template<einhard::LogLevel VERBOSITY, T1> \
inline const einhard::OutputFormatter<VERBOSITY>& operator<<(const einhard::OutputFormatter<VERBOSITY> &out, arg)
#define EINHARD_STREAM_T2(T1, T2, arg1, arg2) \
template<einhard::LogLevel VERBOSITY, T1, T2> \
inline const einhard::OutputFormatter<VERBOSITY>& operator<<(const einhard::OutputFormatter<VERBOSITY> &out, arg1, arg2)

// vim: ts=4 sw=4 tw=100 noet
