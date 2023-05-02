latex_dict = {
    r'\sin': r'\textrm{s}',
    r'\cos': r'\textrm{c}',
    r'\tan': r'\textrm{t}',
    r'{\left(t \right)}': '',
    r'\frac{d}{d t}': r'\dot',
    r'\frac{d^{2}}{d t^{2}}': r'\ddot'
}

def replace_values_in_string(text, args_dict=latex_dict):
    for key in args_dict.keys():
        text = text.replace(key, str(args_dict[key]))
    return text

