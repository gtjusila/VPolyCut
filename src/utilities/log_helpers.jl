using LoggingExtras

function get_logfolder_path()
    logfolder_path = joinpath(dirname(Base.active_project()), "logs")
    if !isdir(logfolder_path)
        mkdir(logfolder_path)
    end
    return logfolder_path
end

function setup_file_logger(filename::String)
    return FormatLogger(open(filename, "w")) do io, args
        # Write the module, level and message only
        println(io, "[", args.level, "] ", args.message)
    end
end