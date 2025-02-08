function prompt_user(;
    message::String="",
    validation::Function=(input) -> true,
    error_message::String="",
    default::String="")
    input = ""
    while (true)
        print(message * " ($(default))" * ": ")
        input = readline()
        input = strip(input)
        if (input == "")
            input = default
        end
        if (validation(input))
            return input
        else
            println(error_message)
        end
    end
end