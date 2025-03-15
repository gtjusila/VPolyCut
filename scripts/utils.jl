function prompt_user(;
    message::String = "",
    validation::Function = (input) -> true,
    parse::Function = (input) -> input,
    error_message::String = "",
    default::String = "")
    input = ""
    while (true)
        print(message * " ($(default))" * ": ")
        input = readline()
        input = strip(input)
        if (input == "")
            input = default
        end
        if (validation(input))
            return parse(input)
        else
            println(error_message)
        end
    end
end