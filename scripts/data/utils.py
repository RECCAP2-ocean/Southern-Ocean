from munch import Munch as _munch


class RECCAP_dict(_munch):
    
    def data(self, dim_name='variable'):
        from xarray import concat, merge, DataArray
        from pandas import Index
        is_array = [isinstance(self[k], DataArray) for k in self]
        is_munch = [isinstance(self[k], self.__class__) for k in self]
        if all(is_array):
            idx = Index(self.keys(), name=dim_name)
            try:
                xda = concat([self[k] for k in self], dim=idx)
            except:
                xda = merge([self[k] for k in self])
            return xda
        else:
            dataarrays = []
            for k in self:
                if isinstance(self[k], self.__class__):
                    dataarrays += self[k].toDataArray(),
            return merge(dataarrays)
        
    def __repr__(self):
        def get_info(obj):
            unit = obj.attrs.get('units', '-')
            dims = ".".join([f"{k}({obj[k].size})" for k in obj.dims])
            unit = f"[UNITS]: {unit}; "
            dims = f"[DIMS]: {dims}"
            info = f"{unit: <32}{dims}"
            return info
        def walk_through_dictionary(d):
            new_dict = {}
            for k, v in d.items():
                if isinstance(v, dict):
                    new_dict[f"{k} {'-'*(76 - len(k))}"] = walk_through_dictionary(v)
                elif isinstance(v, xr.DataArray):
                    new_dict[f"{k: <10}"] = get_info(v)
            return new_dict
        import json
        import xarray as xr
        import re
        import os
        
        printable = walk_through_dictionary(self)
        pretty = f"{str(self.__class__)}"
        pretty += json.dumps(printable, indent=2, sort_keys=True)
        pretty = re.sub('["{},:;]', '', pretty)
        
        m = self.get('not_matched', [])
        if m != []:
            s = "\n    "
            not_matched = s + s.join([os.path.relpath(f) for f in m])
            pretty += f"  not_matched {'='*69}{not_matched}"

        return pretty

        
def read_reccap2_products(flist, fname_specs=[]):
    def flatten_list(list_of_lists):
        if len(list_of_lists) == 0:
            return list_of_lists
        if isinstance(list_of_lists[0], list):
            return flatten_list(list_of_lists[0]) + flatten_list(list_of_lists[1:])
        return list_of_lists[:1] + flatten_list(list_of_lists[1:])
    
    def build_tree(tree_list):
        if tree_list:
            if len(tree_list) > 2:
                return {tree_list[0]: build_tree(tree_list[1:])}
            else:
                xds = xr.open_mfdataset(tree_list[1], 
                                        decode_times=False, 
                                        preprocess=reccap2_dataprep(decode_times=True))
                if len(xds.data_vars) == 1:
                    xds = xds[list(xds.data_vars.keys())[0]]
                return {tree_list[0]: xds}
        return {}
    
    import re
    import xarray as xr
    import mergedeep
    from munch import Munch
    from .preprocess import reccap2_dataprep
    
    flist = flatten_list(flist)

    output = []
    failed = []
    for f in flist:
        if f.endswith('.nc'):
            matches = []
            for spec in fname_specs:
                matches += [m for m in spec if re.findall(spec[m] if spec[m] else m, f)]
            if len(matches) == len(fname_specs):
                tree = [re.sub("[^0-9a-zA-Z]+", "_", m) for m in matches] + [f]
                output += build_tree(tree),
            else:
                failed += f,
    merged = mergedeep.merge({}, *output)            
    merged['not_matched'] = failed
    obj = RECCAP_dict.fromDict(merged)
    
    return obj


