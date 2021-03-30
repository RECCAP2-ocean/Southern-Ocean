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
    merged['failed'] = failed
    obj = RECCAP_dict.fromDict(merged)
    
    return obj


